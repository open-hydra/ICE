# Using ICE

## The case folder

ICE runs in the case folder and reads every path relative to it. Its own files carry
the prefix `part-`.

```
my_case/
├── input.ini                       # the settings of the run
├── INPUT/
│   ├── part-ic.tec                 # the mesh and the initial state
│   ├── part-bc.txt                 # the boundary conditions (ATLAS BCB)
│   ├── part-phase.txt              # the materials (ATLAS); optional with one material
│   ├── part-properties.dat         # the material properties against temperature; optional
│   └── gas.tec                     # the gas field; optional, and it switches the coupling on
└── OUTPUT/                         # create it before the run
    ├── part-field.tec              # the solution, also read by a restart
    ├── part-residual-history.dat   # the residuals
    └── <probe>.txt                 # one file per probe
```

- **The mesh is inside `part-ic.tec`**: its first three variables are the node
  coordinates. There is no separate mesh file. (The `MESH/` folder of the shipped
  cases is only where the mesh generator wrote its files.)
- **The gas file decides the coupling.** If `INPUT/gas.tec` exists, ICE reads the gas
  field and exchanges momentum, heat and, if asked, radiation and mass with it; it
  never changes the gas (one-way coupling). If it does not exist, the cloud moves with
  no exchange at all. `gas-path` changes the folder ICE looks in.
- **`OUTPUT/` must exist** before the run. The `ICE.sh` script of the shipped cases
  creates it.

The formats are described in [Initial Conditions](initial-conditions.md) and
[Boundary Conditions](boundary-conditions.md); the settings in [Input File](input.md).

## Running

```bash
cd my_case
mkdir -p OUTPUT
/path/to/bin/ICE
```

With OpenMP threads:

```bash
export OMP_NUM_THREADS=8
/path/to/bin/ICE
```

With MPI ranks (a build with `--use-mpi`), with or without threads inside each rank:

```bash
export OMP_NUM_THREADS=4
mpirun -np 2 /path/to/bin/ICE
```

Each shipped case has an `ICE.sh` script that does all of this:
`./ICE.sh -m 2 -p 4 solve` runs two ranks of four threads. It also copies `bin/ICE`
into the case when the build is newer, so the case always runs the current build.

[Parallel Execution](parallel.md) explains how to choose ranks and threads.

## The startup report

Before the first iteration, ICE prints what it is about to do. For a coupled run with
evaporation:

```
 Global time step 1-way coupled simulation

 ICE phase model:
 - IG particle families  -->    1

 ICE numerical scheme:
 - Space   --> MUSCL with vanleer flux limiter
 - Time    --> Explicit RK2
 - Drag    --> Schiller-Naumann
 - Heat    --> Ranz-Marshall (Nu = 2 + 0.6 Re^0.5 Pr^1/3)
 - Evap    --> material 1 (A): CEM
 - Interf  --> VLE

 Boundary Conditions:
   Symmetry                       80
   Inflow                         20
   Extrapolation                  300
```

Check it before a long run:

- the first line: a global time step means time accurate, a local one steady state;
  `0-way` means no gas field was found;
- the number of families of each closure is the one you meant;
- `Space` says `MUSCL` if you asked for it, not `first order`;
- every boundary type you expect is listed. A type that is missing has no face
  carrying it.

The exchange models appear only in a coupled run, and the drop's conduction and
boiling lines only when they are switched on.

With MPI, ICE also prints how it split the blocks among the ranks:

```
  MPI partition: 4 blocks over 2 ranks, balance 100.0% of ideal
```

## The solution files

`ic-format` sets the format of the initial state, `sol-format` that of the solution.
They are independent: a case can start from Tecplot ASCII and write VTK.

| Value | Solution file |
|---|---|
| `tecplot ascii` | `OUTPUT/part-field.tec` |
| `tecplot binary` | `OUTPUT/part-field.szplt`; needs a build with TecIO (`-DUSE_TECIO=ON`) |
| `vtk ascii`, `vtk binary`, `vtk raw` | `OUTPUT/part-field.vtm`, with one `.vts` file per block in `OUTPUT/vtk/` |

The file holds the node coordinates, then the variables of every family at the cell
centres: `rho_p1`, `u_p1`, …, `n_p1`, then `rho_p2`, and so on (the list of each
closure is in [Initial Conditions](initial-conditions.md#variable-order)). The ASCII
formats write 17 significant digits, so a value read back is the value written.

`sol-diter` and `sol-dtime` set how often the solution is written. By default each
write replaces the previous file; with `sol-overwrite = false` the files are numbered
instead.

The file also records the simulation time and the counters a restart needs. In
steady-state mode the time is replaced by the iteration count, as a negative number.

## The residual history

`OUTPUT/part-residual-history.dat` gets a row every `res-diter` iterations, and one at
the last iteration. A row holds:

1. the global iteration;
2. the time;
3. to 7. the residuals of the density, the three velocity components and the
   temperature: the L2 norm, over the cells, of their change in one iteration.

With several families, five more columns follow for each family alone, so a row has
$2 + 5\,(1 + n_f)$ values. `res-threshold` reads column 3, the density residual of the
whole domain ([Stopping conditions](input.md#stopping-conditions)).

`shell-diter` sets how often the progress line is printed on screen.

## Restarting

Set `newrun = false` in `[ICE-Parameters]`. ICE then:

- reads its initial state from the solution, `OUTPUT/part-field.*`, in the format of
  `sol-format` (not `ic-format`);
- goes on from the time of that solution;
- appends to the residual history and to the probe files.

A restart of the same kind as the run that wrote the solution — steady after steady,
or time accurate after time accurate — goes on exactly as the uninterrupted run would:
the same iteration count, the same output times and numbers, the same residual
reference. ICE stores these counters in every solution it writes.

Otherwise the restart starts its counters afresh, as a new run does:

- a time-accurate restart from a steady solution starts at $t = 0$;
- a restart from a solution without the counters (written by another program or an
  older ICE) measures `res-threshold` from its own first residual, and says so.

ICE stops instead of guessing in two cases: a Tecplot solution with no
`SOLUTIONTIME` in its zone header (use it as an initial condition, with
`newrun = true`), and a counter that is not finite or out of range. A VTK solution
with no `TIME` is read as time 0.

With grid levels, the restart resumes on the level of the solution. `part-field` is
written on the finest level only, so a restart from it skips the coarse levels.

??? note "Rows written after the solution"
    A run may go on after its last solution, and write residual and probe rows the
    restart will write again. The restart removes them first:

    - from a solution with the counters, any format: the residual rows and the probe
      iteration rows after its iteration count and, time accurate, the probe time rows
      after its time;
    - from a `tecplot ascii` solution without the counters: in a steady restart the rows
      of later iterations; in a time-accurate one the residual rows and the probe time
      rows of later times (the probe iteration rows stay, since such a restart counts
      its iterations from 0).

    In a file written by several runs whose count started again, only the rows after
    the last row not later than the solution go. A last row cut in the middle goes too,
    and a missing final newline is added. From the other solutions nothing is removed:
    delete those rows by hand, and end the file with a newline.

**Numbered solutions.** With `sol-overwrite = false`, rename the solution you want to
restart from to `part-field`. If no `part-field` exists (the run stopped before its
last write), the restart reads the last numbered one.

## Probes

A probe writes the history of a few variables at one point to `OUTPUT/<name>.txt`.
Name it in `[ICE-Probes]`, then give it a section of its own:

```ini
[ICE-Probes]
probe1 = centreline

[centreline]
variables = rho_p1 u_p1 T_p1
position  = 0.5 0.0 0.0
dtime     = 1e-4
```

| Key | Meaning |
|---|---|
| `position` | Physical coordinates: the probe takes the cell whose centre is nearest, over every block |
| `index-position` | `b i j k` instead: a block and the cell indices; a cell outside the mesh stops the run |
| `diter` | A row every that many global iterations, on the finest grid level |
| `dtime` | A row at the first step at or past each multiple of `dtime`, with the time in the first column |
| `variables` | The names of the solution file, as many as you like |

The variables of family `n` are those of the solution files: `rho_p<n>`, `u_p<n>`,
`v_p<n>`, `w_p<n>`, `T_p<n>`, `n_p<n>`, the pressure `P_p<n>` of an IG family, the
tensor `P11_p<n>` … `P33_p<n>` of an AG family, and `f_p<n>`, `chi_p<n>` of a
solidifying family. The same names without the `p` (`rho_<n>`, `P11_<n>`) work too.
`Ts_p<n>` gives the surface temperature of the drops, which differs from `T_p<n>` only
with `liquid-conduction = P2T` or while `boiling = ZGR` boils. Any other name stops
the run at setup.

A new run replaces the probe files. A restart appends to them, and creates those of
probes added since.

## Changing settings during a run

With `ini-diter` set, ICE reads `input.ini` again every that many iterations. The
stopping thresholds and the output frequencies then take effect at once; what is set
up once at the start — the mesh, the scheme, the materials and their properties —
does not change.

A file with an error (a value that breaks its rule, a changed `levels` or `closure`,
a file saved half-written) is ignored with a warning, and the run keeps its previous
settings. With `ini-error = stop`, the run stops instead. See
[Reading again during the run](input.md#reading-again-during-the-run).
