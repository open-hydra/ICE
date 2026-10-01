# ICE strong-scaling harness for monolith

The scripts that produced the numbers of `docs/user/parallel.md`, taken over from
Marco Grossi's campaign of 2026-09-25..28 (`/data10/grossi/tmp/scaling-ice`, jobs
302763, 302933, 302959, 302995; his `results/{FINAL,OPTIMISATION,AUDIT}.md` are the
record of that campaign and stay in the copy at
`/data10/passarani/Desktop/Software/ice-perf/results`). Authorship of the harness,
the case generator (`gencase.f90`, in the copy) and the protocol is his.

## The case

192^3 = 7,077,888 cells in one block, one IG family, particles at rest at
rho_p 1e-2 with P from a 10 m/s dispersion velocity, one-way coupled to a frozen
Taylor-Green carrier (U0 10 m/s), symmetry (300) on all six faces, Tecplot
`.szplt` in and out, 40 iterations, RK2, CFL 0.5, global time stepping (local
stepping diverges by iteration 7 at this size), MUSCL/vanleer, no IRS, no shock
detector, timers on. `INPUT` is the undecomposed mesh; `INPUT-R{2,4,8,12,16,24}`
are MDB cuts with exactly R blocks at 100 % cell balance and 2.1 / 4.2 / 6.2 /
10.4 / 12.5 / 18.8 % ghost cells (`gen-decomp.py`; R = 20 has no cut under its
99 % gate). The decomposition belongs to the RANK count, never the core count.

Local copy: `/data10/passarani/Desktop/Software/ice-perf/case/build` (5.4 GB,
rsynced from Marco's directory on 2026-09-30). On monolith: `/data_fast/gpassarani/
ice-perf/case` (`/scratch` is full).

## What changed on 2026-09-30

* every path and account re-pointed (`ICE_PERF_ROOT`, default
  `/data_fast/gpassarani/ice-perf`; ssh host `monolith`; partition `free`);
* `ice_run` takes EVERY IOK field from the best repetition (the 2026-09 harness
  took the minimum wall/iter but the cycles, phases and clock of the last one);
* `OMP_PROC_BIND=close` set beside `OMP_PLACES=cores` (never set before);
* the job log names the ICE commit, the submodule pins and the dirty-file count
  of the tree that was built (`ICE-commit.txt`, written by `sync-src-monolith.sh`);
* `bench.ini` pins `dt-max` and `tau-factor` (the key did not exist for the
  campaign binaries; `0` reproduces them; neither binds on this case);
* `ice-final.slurm` gains the `blk-*` legs: the R24 cut on one rank and on four,
  the campaign never ran more than one block per rank;
* the `ICE Detail` line of the extended timers is carried into the tsv
  (`d_<region>` columns) when the binary prints it;
* `analyze.py` scores unknown tags in a trailing group instead of aborting, and
  reads ghost percentages from `decompositions.tsv` when given;
* `fetch-results.sh` brings the per-run ICE logs back with the `.out` file (the
  campaign's own run logs exist only under an account we cannot read).

The protocol itself is unchanged: last timer window, one curve per job per node
set, efficiency against the job's own 1-core anchor, minimum over repetitions,
core-cycles beside seconds, no ballast.

## Running a campaign

```sh
# once: the case (5 GB) and, after every change, the tree
test/bench/monolith/deploy-case-monolith.sh
test/bench/monolith/sync-src-monolith.sh hyb          # rsync + build on monolith -> bin/ICE-mpi

# on monolith (bash -l), from ICE_PERF_ROOT
sbatch scripts/ice-smoke.slurm                        # pre-flight, 1 node, ~20 min
sbatch scripts/ice-final.slurm                        # the matrix, 3 nodes, ~12 h
sbatch scripts/ice-omp.slurm                          # placement sweep, 1 node, ~3 h

# back here
test/bench/monolith/fetch-results.sh <jobid>          # -> ice-perf/runs/<jobid>/
python3 test/bench/monolith/collect.py ice-perf/runs/<jobid>/ice-final-<jobid>.out > final.tsv
python3 test/bench/monolith/analyze.py ice-perf/runs/<jobid>/ice-final-<jobid>.out scaling.tsv \
        ice-perf/case/build/decompositions.tsv
```

The reference numbers to compare with are job 302995's (`parallel.md` §6):
14.262 s and 5.469e10 core-cycles per iteration on one core, 80.0 % core-cycle
efficiency at 4x20, 75.3 % at 12x20 on three nodes; noise floor 0.34 % on wall
time and 1.3 % on cycles. A job built from another commit is a new baseline, not
a reproduction: the commit is in its log.

## Traps

* Never rebuild `bin/ICE-mpi` while a job runs: `ice_setup` snapshots the binary
  into `bin.$SLURM_JOB_ID/ICE-mpi` for that reason (job 302988 was lost to a
  mid-flight rebuild).
* The RSS sampler and the occupancy check match on the command name `ICE-mpi`;
  a snapshot under another name records nothing (job 302995 has no memory column).
* `numactl` around `mpirun` does not pin under Slurm; `ICE_WRAP` wraps the binary.
* `analyze.py` needs the `env-C1-1x1` anchor in the same job.
