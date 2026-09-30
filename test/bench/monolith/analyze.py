#!/usr/bin/env python3
"""Score an ICE scaling job.

    python3 analyze.py <job.out> <dest.tsv> [<decompositions.tsv>]

Two efficiencies are reported for every point, and they answer different
questions:

  eff_wall  = (T(1)/N) / T(N)     what the user waits for
  eff_cyc   = C(1) / C(N)         what the code does

C is total core-cycles per iteration summed over every rank and thread. Perfect
parallelism keeps it constant, so eff_cyc has the chip's frequency scaling
divided out -- one active core turbos to 3.83 GHz on this node while eighty sit
near 2.67, and that factor is the machine, not ICE. eff_wall is always the
smaller of the two; the ratio between them is exactly the clock ratio, which is
a useful check that the instrumentation is measuring what it claims.

eff_cyc_gn additionally credits the decomposition's ghost cells, which are real
extra work MDB adds and which the 1-block anchor never does. It is the fairest
statement of how well the parallelisation itself holds up, and it is the number
to quote when comparing ICE against another code on another mesh. The ghost
percentages come from decompositions.tsv when given (the 2026-09 campaign's
numbers are the fallback).

Tags outside the known groups are still scored, in a trailing "other" group,
so a job with extra legs does not abort the table.
"""
import re, sys, collections

GHOST_DEFAULT = {1: 0.0, 2: 2.1, 4: 4.2, 8: 6.2, 12: 10.4, 16: 12.5, 24: 18.8}

GROUPS = [("envelope, hybrid 4 ranks/node",
           ["env-C1-1x1", "env-C4-4x1", "env-C8-4x2", "env-C16-4x4",
            "env-C20-4x5", "env-C40-4x10", "env-C80-4x20",
            "env-C160-8x20", "env-C240-12x20"]),
          ("pure MPI", ["mpi-C2-2x1", "mpi-C4-4x1", "mpi-C8-8x1",
                        "mpi-C16-16x1", "mpi-C24-24x1"]),
          ("pure OpenMP", ["omp-C2-1x2", "omp-C4-1x4", "omp-C8-1x8",
                           "omp-C16-1x16", "omp-C20-1x20", "omp-C40-1x40",
                           "omp-C80-1x80"]),
          ("layout at 80 cores", ["lay-C80-4x20", "lay-C80-8x10", "lay-C80-16x5"]),
          ("socket coverage at 20 cores",
           ["sock-C20-1x20-packed", "sock-C20-1x20-spread", "sock-C20-4x5-spread"]),
          ("many blocks per rank (R24 cut)", ["blk-C80-1x80-R24", "blk-C80-4x20-R24"]),
          ("anchor", ["anchor-4x20"])]


def load(path):
    rows = []
    for line in open(path):
        if not line.startswith("IOK "):
            continue
        d = dict(re.findall(r"(\w+)=(\S*)", line))
        d["cores"] = int(d["ranks"]) * int(d["threads"])
        for k in ("wall_iter", "cycles_iter", "ghz", "rss_rank_gb",
                  "rss_node_gb", "setup_s", "halo_pct", "collective_wait_pct",
                  "exchange_wait_pct", "compute_spread_pct"):
            d[k] = float(d[k]) if d.get(k) else None
        rows.append(d)
    return rows


def load_ghost(path):
    """set -> ghost_pct from gen-decomp.py's decompositions.tsv."""
    ghost = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("set"):
            continue
        f = line.split()
        if len(f) >= 6:
            ghost[f[0]] = float(f[5])
    return ghost


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    rows = load(sys.argv[1])
    by = {r["tag"]: r for r in rows}
    anchor = by["env-C1-1x1"]
    T1, C1 = anchor["wall_iter"], anchor["cycles_iter"]
    ghost_by_set = load_ghost(sys.argv[3]) if len(sys.argv) > 3 else {}

    def ghost_pct(r):
        if r.get("set") in ghost_by_set:
            return ghost_by_set[r["set"]]
        return GHOST_DEFAULT.get(int(r["ranks"]), 0.0)

    def score(r):
        N = r["cores"]
        gh = 1.0 + ghost_pct(r) / 100.0
        return (T1 / N / r["wall_iter"] * 100,
                C1 / r["cycles_iter"] * 100,
                C1 * gh / r["cycles_iter"] * 100,
                r["ghz"] / N)

    dest = sys.argv[2]
    out = open(dest, "w")
    # The job id comes from the log, not from a literal: a hard-coded job number
    # in the header of every table is how a post-fix curve ends up labelled with
    # the baseline's job.
    job = node = commit = binary = "unknown"
    for line in open(sys.argv[1]):
        if line.startswith("JOB="):
            m = re.search(r"id=(\S+)", line); job = m.group(1) if m else job
            m = re.search(r"nodes=(\S+)", line); node = m.group(1) if m else node
        elif line.startswith("BINARY"):
            m = re.search(r"md5=(\S+)", line); binary = m.group(1) if m else binary
        elif line.startswith("COMMIT"):
            commit = line.split(None, 1)[1].strip()
    out.write("# ICE scaling, job %s, nodes %s, binary %s\n" % (job, node, binary))
    out.write("# commit %s\n" % commit)
    out.write("# case 192^3 = 7,077,888 cells, IG closure, 40 iterations\n")
    out.write("# anchor: 1 core, %.3f s/iter, %.4e core-cycles/iter\n" % (T1, C1))
    out.write("# eff_cyc_gn credits the decomposition's ghost cells; see analyze.py\n")
    out.write("tag\tcores\tranks\tthreads\twall_iter_s\teff_wall_pct\t"
              "cycles_iter\teff_cyc_pct\teff_cyc_gn_pct\tclock_ghz_per_core\t"
              "halo_pct\tcollective_wait_pct\tcompute_spread_pct\trss_node_gb\tsetup_s\n")

    seen = set()
    groups = list(GROUPS)
    other = [t for t in by if not any(t in tags for _, tags in GROUPS)]
    if other:
        groups.append(("other", sorted(other)))

    for title, tags in groups:
        tags = [t for t in tags if t in by]
        if not tags:
            continue
        print("\n== %s ==" % title)
        print("%-22s %5s %8s %7s %11s %7s %7s %6s" %
              ("tag", "cores", "s/iter", "eff_w%", "cycles/iter", "eff_c%", "eff_gn%", "GHz/c"))
        for t in tags:
            r = by[t]
            seen.add(t)
            ew, ec, eg, gh = score(r)
            print("%-22s %5d %8.4f %7.1f %11.4e %7.1f %7.1f %6.2f" %
                  (t, r["cores"], r["wall_iter"], ew, r["cycles_iter"], ec, eg, gh))

            # Any of the sampled columns can be absent -- the RSS sampler catches
            # nothing on a point that finishes inside its 5 s poll, and a
            # single-rank run reports no spread. Write an empty field rather
            # than letting one missing number abort the whole campaign's tsv.
            def fmt(v, spec):
                return "" if v is None else spec % v
            out.write("\t".join([
                t, "%d" % r["cores"], str(r["ranks"]), str(r["threads"]),
                "%.4f" % r["wall_iter"], "%.1f" % ew,
                "%.4e" % r["cycles_iter"], "%.1f" % ec, "%.1f" % eg, "%.3f" % gh,
                fmt(r["halo_pct"], "%.1f"),
                fmt(r["collective_wait_pct"], "%.1f"),
                fmt(r["compute_spread_pct"], "%.1f"),
                fmt(r["rss_node_gb"], "%.2f"),
                fmt(r["setup_s"], "%.1f")]) + "\n")
    out.close()

    # A campaign whose RSS sampler caught nothing must SAY so rather than divide
    # by zero or print a column of zeros that reads like a measurement.
    mem = ["env-C1-1x1", "mpi-C2-2x1", "mpi-C4-4x1", "mpi-C8-8x1",
           "mpi-C16-16x1", "mpi-C24-24x1"]
    if all(t in by and by[t]["rss_node_gb"] for t in mem):
        print("\n== memory: every rank holds the whole domain ==")
        for t in mem:
            r = by[t]
            n = int(r["ranks"])
            print("  %2d ranks: %6.2f GB/node  %.2f GB/rank" %
                  (n, r["rss_node_gb"], r["rss_node_gb"] / n))
        print("  187 GB/node => cap ~%d ranks/node at this size" %
              int(187 / (by["mpi-C24-24x1"]["rss_node_gb"] / 24)))
    else:
        print("\n== memory: NOT MEASURED in this campaign (rss sampler empty) ==")
    print("\nwrote %s" % dest)


main()
