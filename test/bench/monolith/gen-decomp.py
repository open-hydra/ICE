#!/usr/bin/env python3
"""Build one decomposition per MPI RANK count for the ICE scaling campaign.

The decomposition belongs to the rank count, never the core count: a hybrid
4x20 run is four ranks and wants the 4-block grid, not an 80-block one.
Indexing the grid by core count is what broke the OpenMP curve in the reference
MOSE campaign, and the hydra-MF campaign inherited the rule from it.

ICE reads two grids that must be cut identically -- the condensed phase
(part-ic.szplt, with part-bc.txt) and the frozen carrier (gas.szplt) -- so both
go through one MDB run as Phase1 and Phase2. They share the mesh, so the same
parameters give the same cut; the script checks that they did.

`min-cells` is the lever that stops MDB slabbing a block into a thin sheet. On
a cube MDB's greedy longest-direction rule mostly does the right thing on its
own, but the awkward rank counts (20 = 2^2 x 5, 24, 12) are worth sweeping, so
each rank count is tried at several values and the lowest ghost overhead that
still reaches the rank count at full balance wins.

Paths come from the environment (2026-09-30): ICE_PERF_CASE_LOCAL is the case
directory holding INPUT/ and mdb.tmpl, ATLASDIR the ATLAS tree whose bin/MDB is
run. Marco Grossi's original lives in /data10/grossi/tmp/scaling-ice/scripts.
"""
import os, re, shutil, subprocess, sys

CASE = os.environ.get("ICE_PERF_CASE_LOCAL",
                      "/data10/passarani/Desktop/Software/ice-perf/case/build")
ATLASDIR = os.environ.get("ATLASDIR", "/data10/passarani/Desktop/Software/hydra/utils/ATLAS")
MDB = os.environ.get("MDB", os.path.join(ATLASDIR, "bin", "MDB"))
RANKS = [2, 4, 8, 12, 16, 20, 24]
MINCELLS = [24, 48, 96]

# The min-cells that won the sweep, per rank count. MDB's cut is geometric --
# it depends on the mesh, not on the field values -- so once the sweep has been
# run on a given mesh its answers stay valid for any regenerated IC on that same
# mesh. `--known` re-cuts with these directly, 7 MDB runs instead of 21.
# R=20 is absent on purpose: 192 = 2^6 x 3 has no factor of 5, and MDB's best
# 20-way split under the 99 % gate reaches only 80% balance, so the campaign
# skips that rank count. (A 5x4x1 cut at ~98.5 % exists with min-cells 33-38;
# Phase A of the parallelization plan is where MDB learns to find it.)
KNOWN = {2: 24, 4: 96, 8: 96, 12: 48, 16: 48, 24: 24}

TMPL = open(os.path.join(CASE, "mdb.tmpl")).read()


def run_mdb(R, MC, out):
    ini = TMPL.replace("@RANKS@", str(R)).replace("@MINCELLS@", str(MC)).replace("@OUT@", out)
    path = os.path.join(CASE, "mdb-R%d.ini" % R)
    open(path, "w").write(ini)
    env = dict(os.environ, ATLASDIR=ATLASDIR)
    p = subprocess.run([MDB, "-i", path], cwd=CASE, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode(errors="replace")


def parse(log):
    """(blocks, balance, ghost) per phase; MDB prints one summary per phase."""
    nb = [int(m) for m in re.findall(r"New blocks\s+(\d+)", log)]
    bal = [float(m) for m in re.findall(r"balance\s+([0-9.]+)%", log)]
    gh = [float(m) for m in re.findall(r"overhead\s+([0-9.]+)%", log)]
    return nb, bal, gh


def main():
    known = "--known" in sys.argv
    rows = []
    for R in (sorted(KNOWN) if known else RANKS):
        best = (KNOWN[R], None, None, None) if known else None
        for MC in ([] if known else MINCELLS):
            out = ".try-R%d" % R
            shutil.rmtree(os.path.join(CASE, out), ignore_errors=True)
            os.makedirs(os.path.join(CASE, out), exist_ok=True)
            rc, log = run_mdb(R, MC, out)
            nb, bal, gh = parse(log)
            ok = rc == 0 and len(nb) == 2 and nb[0] == nb[1] == R and min(bal) >= 99.0
            print("    R=%-3d min-cells=%-3d -> blocks=%s balance=%s ghost=%s %s"
                  % (R, MC, nb, bal, gh, "" if ok else "(rejected)"), flush=True)
            if ok and (best is None or gh[0] < best[1]):
                best = (MC, gh[0], nb[0], bal[0])
            shutil.rmtree(os.path.join(CASE, out), ignore_errors=True)
        if best is None:
            print("  R=%d: NO ADMISSIBLE DECOMPOSITION" % R, flush=True)
            continue
        MC, gh, nb, bal = best
        out = "INPUT-R%d" % R
        shutil.rmtree(os.path.join(CASE, out), ignore_errors=True)
        os.makedirs(os.path.join(CASE, out), exist_ok=True)
        rc, log = run_mdb(R, MC, out)
        if known:
            nbk, balk, ghk = parse(log)
            if rc != 0 or len(nbk) != 2 or nbk[0] != R or min(balk) < 99.0:
                print("  R=%-3d REJECTED on re-cut: rc=%d blocks=%s balance=%s"
                      % (R, rc, nbk, balk), flush=True)
                continue
            nb, bal, gh = nbk[0], balk[0], ghk[0]
        shutil.copy(os.path.join(CASE, "INPUT/part-phase.txt"),
                    os.path.join(CASE, out, "part-phase.txt"))
        # ICE never reads the gas BC file; MDB only needs it to cut Phase2.
        for junk in ("gasbc.txt",):
            j = os.path.join(CASE, out, junk)
            if os.path.exists(j):
                os.remove(j)
        print("  R=%-3d -> %s : %d blocks, min-cells=%d, ghost %.1f%%, balance %.1f%%"
              % (R, out, nb, MC, gh, bal), flush=True)
        rows.append((R, nb, MC, bal, gh))

    dest = os.path.join(CASE, "decompositions.tsv")
    with open(dest, "w") as f:
        f.write("# Decomposition used at each MPI rank count, from gen-decomp.py.\n")
        f.write("#\n")
        f.write("# The decomposition belongs to the RANK count, not the core count: a hybrid\n")
        f.write("# 4x20 run is four ranks and uses R4.\n")
        f.write("#\n")
        f.write("# ghost_pct is the extra cell count the decomposition adds over the\n")
        f.write("# undecomposed mesh. It is real work, so a strong-scaling curve taken across\n")
        f.write("# these sets mixes parallel inefficiency with decomposition cost; that is why\n")
        f.write("# the analysis reports both a raw and a ghost-normalised efficiency.\n")
        f.write("set\tranks\tblocks\tmin_cells\tbalance_pct\tghost_pct\n")
        f.write("INPUT\t1\t1\t-\t100.0\t0.0\n")
        for R, nb, MC, bal, gh in rows:
            f.write("INPUT-R%d\t%d\t%d\t%d\t%.1f\t%.1f\n" % (R, R, nb, MC, bal, gh))
    print("wrote %s" % dest)


if __name__ == "__main__":
    main()
