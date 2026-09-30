#!/usr/bin/env python3
"""Turn an ICE scaling job's .out file into a tsv, one row per measured point.

Reads the IOK lines the harness emits (one per point, every field from the
best repetition, scored on the LAST timer window) plus the IOCC occupancy
lines, and writes a tab-separated table with a provenance header.

    python3 collect.py <job.out> [<job.out> ...] > runs/<job>/final.tsv

The optional `detail=` field (extended timers, 2026-09-30) is split into one
column per region, named `d_<region>`; older logs leave those columns empty.
"""
import os, re, sys

FIELDS = ["tag", "set", "ranks", "threads", "ppn", "reps", "rep", "wall_iter",
          "rank_min_s", "rank_avg_s", "imbalance_pct",
          "exchange_wait_pct", "collective_wait_pct",
          "source_s", "source_pct", "flux_s", "flux_pct", "halo_s", "halo_pct",
          "compute_max_s", "compute_min_s", "compute_mean_s", "compute_spread_pct",
          "cycles_iter", "ghz", "blocks", "balance",
          "rss_rank_gb", "rss_node_gb", "setup_s", "loop_s", "iters"]

HEADER = """\
# ICE strong-scaling campaign, generated {date}
# case: 7,077,888 cells (192^3), IG closure, one-way coupled to a frozen
#       Taylor-Green carrier; 40 iterations, RK2, MUSCL/vanleer
# metric: wall_iter = seconds per iteration, LAST timer window, min over reps;
#   every other column comes from that same best repetition (rep)
# cycles_iter = total core-cycles per iteration, summed over every rank and
#   thread. Parallel efficiency from seconds also charges the solver for the
#   chip's frequency scaling: one active core turbos to ~3.8 GHz here while
#   eighty sit near ~2.6 GHz, a factor that has nothing to do with the code.
#   Cycles are immune -- perfect parallelism keeps total core-cycles per
#   iteration constant whatever frequency each core runs at, so C(1)/C(N) is
#   efficiency with the clock already divided out. Prefer it over seconds.
# 'set' is the decomposition, which belongs to the RANK count, not the core
#   count: a hybrid 4x20 run is four ranks and uses INPUT-R4.
# rss_rank_gb is the largest single rank. ICE gives every rank the whole domain
#   (Wrap_Setup: "every rank keeps the whole domain but updates only its own
#   blocks"), so this does NOT fall as ranks are added, and the node's memory
#   is what caps pure MPI -- not its core count.
# setup_s is wall time outside the 40 iterations: grid read, decomposition and
#   the final write. Every rank reads the whole grid, so it grows with ranks.
# clock: no ballast. Seconds are as-measured on an otherwise idle node, which
#   is what a user gets; cycles_iter is the clock-independent statement.
# imbalance_pct = (slowest rank - mean rank)/mean, over the iteration time.
# exchange_wait_pct / collective_wait_pct = share of the summed rank time spent
#   blocked on the halo exchange and on collectives respectively.
# compute_spread_pct = the same spread over compute alone (iteration minus both
#   waits), which separates "this rank had more work" from "this rank waited".
# d_<region> = seconds per iteration of that timer region (max over ranks),
#   from the ICE Detail line; empty on binaries without the extended timers.
{prov}#
"""


def split_detail(s):
    """'dt;1.2E-3;copy;3.4E-4;...' -> {'dt': '1.2E-3', ...} (tolerant of junk)."""
    parts = [p for p in s.split(";") if p]
    out = {}
    for i in range(0, len(parts) - 1, 2):
        key = re.sub(r"[^A-Za-z0-9_]", "", parts[i]).lower()
        if key:
            out[key] = parts[i + 1]
    return out


def main():
    rows, prov = [], []
    for path in sys.argv[1:]:
        job = node = binary = commit = ""
        for line in open(path, errors="replace"):
            if line.startswith("JOB="):
                m = dict(re.findall(r"(\w+)=(\S+)", line))
                job, node = m.get("JOB", ""), m.get("nodes", "")
            elif line.startswith("BINARY"):
                binary = re.search(r"md5=(\S+)", line).group(1)
            elif line.startswith("COMMIT"):
                commit = line.split(None, 1)[1].strip()
            elif line.startswith("IOK "):
                d = dict(re.findall(r"(\w+)=(\S+)", line))
                d["_src"] = os.path.basename(path)
                d["_node"] = node
                if "detail" in d:
                    for k, v in split_detail(d.pop("detail")).items():
                        d["d_" + k] = v
                rows.append(d)
        prov.append("# %s: job=%s node=%s binary=%s commit=%s\n"
                    % (os.path.basename(path), job, node, binary, commit))

    detail_cols = sorted({k for d in rows for k in d if k.startswith("d_")})
    cols = FIELDS + detail_cols

    import datetime
    sys.stdout.write(HEADER.format(date=datetime.date.today().isoformat(),
                                   prov="".join(prov)))
    sys.stdout.write("node\t" + "\t".join(cols) + "\n")
    for d in rows:
        sys.stdout.write(d.get("_node", "") + "\t"
                         + "\t".join(d.get(f, "") for f in cols) + "\n")
    sys.stderr.write("collected %d points from %d file(s)\n" % (len(rows), len(sys.argv) - 1))


if __name__ == "__main__":
    main()
