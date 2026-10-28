# Using ICE

## Case directory layout

ICE runs in the case directory and resolves every path relative to it. All of its
files are prefixed `part-`.

```
my_case/
├── input.ini                       # the only configuration file
├── INPUT/
│   ├── part-ic.tec                 # initial condition, and the mesh with it
│   ├── part-bc.txt                 # boundary condition table (ATLAS BCB)
│   ├── part-properties.dat         # optional rho(T), cs(T) table
│   └── gas.tec                     # optional gas field; its presence enables coupling
└── OUTPUT/                         # created by ICE
    ├── part-field.tec              # solution, and the restart file
    ├── part-residual-history.dat   # residual per iteration
    └── <probe>.txt                 # one per probe, if any
```

`OUTPUT/` must exist before the run; the `ICE.sh` script in each test case creates it.
The mesh is carried in `part-ic.tec` itself — the node coordinates are the first three
variables of the file — so there is no separate mesh input. A `MESH/` directory appears
in the shipped cases only because the mesh generator wrote it there.

Whether the run is coupled is decided by one thing: whether `INPUT/gas.tec` exists. If
it does, drag, heat and, with a non-zero `emissivity`, radiation are active and the gas
field is read but never modified (1-way coupling); if it does not, the cloud is transported
with no source terms at all (0-way). `gas-path` moves the directory ICE looks in.

## Workflow

1. **Prepare the inputs** — mesh and initial condition in `INPUT/part-ic.tec`
   ([format](initial-conditions.md)), the boundary table in `INPUT/part-bc.txt`
   ([format](boundary-conditions.md)), optionally a property table and a gas field.
2. **Write `input.ini`** — at minimum one `[ICE-FamilyN]` section naming the closure,
   and a stopping condition. See the [input file](input.md) page and the full
   [parameter reference](registry.md).
3. **Run** — `bin/ICE` from the case directory.
4. **Post-process** — open `OUTPUT/part-field.tec` (or the VTK equivalent) in Tecplot,
   ParaView or VisIt.

## Running

```bash
cd my_case
mkdir -p OUTPUT
/path/to/bin/ICE
```

With OpenMP:

```bash
export OMP_NUM_THREADS=8
/path/to/bin/ICE
```

With MPI (the build must use `--use-mpi`), optionally with OpenMP threads inside each
rank:

```bash
export OMP_NUM_THREADS=4
mpirun -np 2 /path/to/bin/ICE
```

The `ICE.sh` script in each test case wraps all three: `./ICE.sh -m 2 -p 4 solve` runs
on two ranks of four threads. It also copies the master `bin/ICE` into the case when
the master is newer, so a case always runs against the current build.

### What the startup report tells you

```
 ICE phase model:
 - AG particle families  -->    1

 ICE numerical scheme:
 - Space   --> MUSCL with MC flux limiter
 - Time    --> Explicit Runge-Kutta 2
 - Drag    --> Stokes
 - Heat    --> Ranz-Marshall (Nu = 2 + 0.6 Re^0.5 Pr^1/3)
 - Evap    --> CEM

 Boundary Conditions:
   Extrapolation                  2
```

The exchange models appear only in a coupled run, and the interface and blowing lines
only when an evaporation model is selected. Worth checking before a long run: the
closure count matches the families you meant to declare, the scheme line says `MUSCL`
rather than `I order` if you asked for it, and every boundary type you expect appears
with the right count. A type that is absent
from the list has no faces carrying it.

### How MPI divides the work

ICE distributes **whole blocks** over the ranks, largest first, each to the rank with
the least work so far. At startup it prints how even the split is:

```
  MPI partition: 4 blocks over 2 ranks, balance 100.0% of ideal
  MPI halo: 400 cells exchanged per ghost fill
  MPI halo: 100 connection faces whose flux crosses ranks
```

- A block is never split, so ranks beyond the number of blocks sit idle. To use more
  ranks, split the mesh into more blocks.
- Every rank reads the whole mesh and solution, and updates only its own blocks. Before
  each ghost-cell fill — that is, once per Runge-Kutta stage — the interior cells that
  a connection (`101`/`201`) or chimera (`102`) face reads from a block owned by
  another rank are sent to it. The flux through a connection face is computed once, by
  the rank of one of its two blocks, and sent to the rank of the other, which first sends
  the shock-detector weight of its cell.
- Rank 0 gathers the blocks and writes the output files.

The result does not depend on the number of ranks: the solution is bit-for-bit the same
as a serial run. See [Multi-block and MPI](../vv/multiblock-mpi.md).

For how to choose the rank/thread split, where to place the threads, and how well ICE
scales on one node and across three, see [Parallel Execution](parallel.md).

## Input and output formats

Two keys, both a writer and a mode: `ic-format` names the format of
`INPUT/part-ic.*`, and `sol-format` the format of `OUTPUT/part-field.*`. They are
independent — a case may be given in Tecplot ASCII and written in VTK.

| Value | File |
|---|---|
| `tecplot ascii` | `part-field.tec` |
| `tecplot binary` | `part-field.szplt`; needs a TecIO-enabled build |
| `vtk ascii` / `vtk binary` / `vtk raw` | `part-field.vtm` plus one `.vts` per block under `vtk/` |

The file holds the node coordinates followed by the primitive variables of every
family, cell-centred, named `rho_p1`, `u_p1`, … `n_p1`, `rho_p2`, …, a solidifying family adding `f_p1` and
`chi_p1` after `n_p1` (a probe takes the same names), every value to 17 significant
digits in the ASCII formats. The solution time is recorded in the zone header
(`SOLUTIONTIME`) and in the `TIME` field data of a vtk file; in steady-state mode it
carries the global iteration count instead, negated. Next to it every solution stores the
counters a restart takes back (see [Restarting](#restarting)): data set auxiliary data
records after the zones of a Tecplot file (`DATASETAUXDATA CYCLE = "40"`, ...), field data
of every block file of a vtk one.

`sol-diter` and `sol-dtime` control how often the solution is written. With
`sol-overwrite = false` each write appends a counter to the name, giving a numbered
series instead of one file that is repeatedly replaced.

`shell-diter` sets how often the progress line is printed and `res-diter` how often a
row is appended to `OUTPUT/part-residual-history.dat`.

## Restarting

Set `newrun = false`. ICE then reads its initial state from `OUTPUT/part-field.*`
instead of `INPUT/part-ic.*`, goes on from the time of that solution, and continues
appending to the residual history. The restart file is a solution ICE wrote itself, so
it is `sol-format` — not `ic-format` — that names it.

The time of the solution is its `SOLUTIONTIME` (Tecplot, ASCII or binary) or its `TIME`
(vtk). A negative one marks a steady solution, and a steady restart goes on from its
iteration count; the Tecplot files of earlier versions, which hold the count positive with
`STRANDID = 0` in their zone headers, are read the same way. A time-accurate restart goes
on from the time of a time-accurate solution, and starts at $t = 0$ from a steady one. A
Tecplot solution whose zone header has no `SOLUTIONTIME` cannot tell where to go on: the
restart stops (use the file as the initial condition instead, with `newrun = true`). A vtk
solution without `TIME` reads as time 0.

Every solution also stores the counters that make a restart go on as the uninterrupted
run would: `CYCLE`, the global iteration count (the first column of the residual
history), `ICE_ITERATION`, the iterations of the grid level of the solution, the finest
for `part-field` (which the CFL ramp of `cfl-rise-threshold` and the `sol-diter` and
`res-diter` cadences count),
`ICE_OUTPUT_NUMBER`, the last number of the series of `sol-overwrite = false`,
`ICE_GRID_LEVEL`, the grid level being solved, `ICE_LAST_OUTPUT_TIME`, the time of the
last `sol-dtime` solution, and `ICE_RES_REFERENCE`, the unit of `res-threshold`. A
restart of the same kind as the run that wrote the solution, steady after steady or
time-accurate after time-accurate, takes them back, and writes the solutions, the
residual history and the probe rows that the uninterrupted run would have written; every
run writes a residual row at its last iteration, so a restart that does not fall on a
multiple of `res-diter` adds that one row. A restart of the other kind, or from a solution
without them (written by an earlier version of ICE, or by another program), starts them
as a new run does: a steady one then measures `res-threshold` from its own first non-zero
residual, and says so. A multigrid run resumes on the grid level of its solution;
`part-field` is written on the finest level only, so a restart from it skips the coarse
levels. A counter that is not finite, or a count out of range, stops the run.

A run that went on past the solution the restart reads (it stopped after its last
output) left in the residual history and in the probe files rows that the restart writes
again, and the restart drops them first. From a solution that holds the counters, in any
format, it drops the rows of the iterations after its `CYCLE` (the residual rows and the
iteration rows of a probe) and, time-accurate, the time rows of a probe after its time.
From a `tecplot ascii` solution without them, a steady restart drops the rows of later
iterations, a time-accurate one the residual rows and the time-based probe rows of later
times (the iteration rows of a probe stay: such a restart counts its iterations from 0).
In a file written by several runs whose count started again, the rows before the last one
not later than the solution stay. A last row cut in the middle goes too, and a file whose
last line has no newline gets one. From the other solutions nothing is dropped: remove
those rows by hand, and end the file with a newline.

With `sol-overwrite = false` the solutions form a numbered series; restarting from a
particular one means renaming it back to `part-field`. A restart that finds no
`part-field` (a run stopped before its last write) reads the last of the series.

## Probes

A probe writes the history of a few variables at one point to `OUTPUT/<name>.txt`.
Declare it in two steps — a reference in `[ICE-Probes]` and a section of its own:

```ini
[ICE-Probes]
probe1 = centreline

[centreline]
variables      = rho_p1 u_p1 T_p1
position       = 0.5 0.0 0.0
dtime          = 1e-4
```

Give either `position` (physical coordinates: the probe takes the cell whose centre is
nearest, over every block) or
`index-position` as `b i j k`; an `index-position` that is not a cell of the mesh stops
the run. `dtime` and `diter` set how often the probe samples: a row every `diter`
iterations of the global count (the first column of the residual history), written on the
finest grid level, and a row at the first step at or past each multiple of `dtime` of
physical time, from the first, with the time in its first column (one row for a step past
several multiples; a step can write both rows). A restart that goes on in time writes the
time rows the uninterrupted run would.

The variables of family `n` are those of the solution files, `rho_p<n>`, `u_p<n>`,
`v_p<n>`, `w_p<n>`, `T_p<n>`, `n_p<n>`, and `f_p<n>`, `chi_p<n>` for a solidifying family;
the same names without the `p` (`rho_<n>`, …) are accepted too. Any other name stops the
run at set-up.

A new run replaces the probe files; a restart appends to them, after the rows written past
its solution are dropped (see [Restarting](#restarting)), and creates any that are
missing (for example a probe added before the restart).

## Re-steering a running simulation

If `ini-diter` is set, ICE re-reads `input.ini` every that many iterations, so
thresholds (the iteration budgets included) and output frequencies can be changed while
the run is in progress. Values that are consumed once at setup — the mesh, the scheme,
the materials and their properties — are not affected. An edit with an error (a value
that breaks its rule, a changed `levels` or `closure`, a file saved half-written) is not
taken by default: ICE warns and the run keeps every value it had; with
`ini-error = stop` it stops instead (see [Re-reading at runtime](input.md#re-reading-at-runtime)).
