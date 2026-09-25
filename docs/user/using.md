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
it does, drag, heat and radiation are active and the gas field is read but never
modified (1-way coupling); if it does not, the cloud is transported with no source
terms at all (0-way). `gas-path` moves the directory ICE looks in.

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
```

- A block is never split, so ranks beyond the number of blocks sit idle. To use more
  ranks, split the mesh into more blocks.
- Every rank reads the whole mesh and solution, and updates only its own blocks. Before
  each ghost-cell fill — that is, once per Runge-Kutta stage — the interior cells that
  a connection (`101`/`201`) or chimera (`102`) face reads from a block owned by
  another rank are sent to it.
- Rank 0 gathers the blocks and writes the output files.

The result does not depend on the number of ranks: the solution is bit-for-bit the same
as a serial run. See [Multi-block and MPI](../vv/multiblock-mpi.md).

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
family, cell-centred, named `rho_p1`, `u_p1`, … `n_p1`, `rho_p2`, … The solution time
is recorded in the zone header; in steady-state mode it carries the iteration count
instead, negated.

`sol-diter` and `sol-dtime` control how often the solution is written. With
`sol-overwrite = false` each write appends a counter to the name, giving a numbered
series instead of one file that is repeatedly replaced.

`shell-diter` sets how often the progress line is printed and `res-diter` how often a
row is appended to `OUTPUT/part-residual-history.dat`.

## Restarting

Set `newrun = false`. ICE then reads its initial state from `OUTPUT/part-field.*`
instead of `INPUT/part-ic.*`, picks up the simulation time from the zone header, and
continues appending to the residual history. The restart file is a solution ICE wrote
itself, so it is `sol-format` — not `ic-format` — that names it.

With `sol-overwrite = false` the solutions form a numbered series; restarting from a
particular one means renaming it back to `part-field`.

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

Give either `position` (physical coordinates, and ICE finds the nearest cell) or
`index-position` as `b i j k`. `dtime` and `diter` set how often the probe samples.

## Re-steering a running simulation

If `ini-diter` is set, ICE re-reads `input.ini` every that many iterations, so
thresholds and output frequencies can be changed while the run is in progress. Values
that are consumed once at setup — the closures, the mesh, the scheme, the materials and
their properties — are not affected.
