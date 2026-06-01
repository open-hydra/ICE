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

With MPI:

```bash
mpirun -np 4 /path/to/bin/ICE
```

With OpenMP:

```bash
export OMP_NUM_THREADS=8
/path/to/bin/ICE
```

## Restarting

ICE writes restart (backup) files at intervals specified by `bck_diter` or `bck_dtime` in `input.ini`. To restart from a backup, set the appropriate restart flag in `input.ini` and re-run.

---
