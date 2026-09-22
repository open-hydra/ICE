# Using ICE

## Workflow

A typical ICE simulation follows this workflow:

1. **Prepare input files** – mesh, initial conditions, boundary conditions, and `input.ini`.
2. **Configure** – edit `input.ini` to set solver parameters, output format, and particle group properties.
3. **Run** – execute `bin/ICE` from the case directory.
4. **Post-process** – open the output files in Tecplot or ParaView.

## Case directory layout

```
my_case/
├── input.ini          # Main configuration file
├── grid/              # Multi-block structured grid files
├── ic/                # Initial condition files (one per block per group)
├── bc/                # Boundary condition tables
└── output/            # Written by ICE at runtime
```

## Running ICE

```bash
cd my_case
/path/to/bin/ICE
```

With OpenMP:

```bash
export OMP_NUM_THREADS=8
/path/to/bin/ICE
```

With MPI (the build must use `--use-mpi`), optionally combined with OpenMP threads in
each rank:

```bash
export OMP_NUM_THREADS=4
mpirun -np 2 /path/to/bin/ICE
```

The `ICE.sh` script in each test case does the same with `./ICE.sh -m 2 -p 4 solve`.

### How MPI divides the work

ICE distributes **whole blocks** over the ranks, largest first, each to the rank with the
least work so far. At startup it prints how even the split is:

```
  MPI partition: 4 blocks over 2 ranks, balance 100.0% of ideal
  MPI halo: 400 cells exchanged per ghost fill
```

- A block is never split, so ranks beyond the number of blocks sit idle. To use more
  ranks, split the mesh into more blocks.
- Every rank reads the whole mesh and solution, and updates only its own blocks. Before
  each ghost-cell fill, the interior cells that a connection (`101`/`201`) or chimera
  (`102`) face reads from a block owned by another rank are sent to it.
- Rank 0 collects the blocks and writes the output files.

The result does not depend on the number of ranks: the solution is bit-for-bit the same as
a serial run.

## Restarting

ICE writes restart (backup) files at intervals specified by `bck_diter` or `bck_dtime` in `input.ini`. To restart from a backup, set the appropriate restart flag in `input.ini` and re-run.

---
