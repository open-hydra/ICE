# Quick Start

This page runs one of the shipped cases, end to end, after a successful
[installation](installation.md).

## 1. Check the build

```bash
ls bin/ICE
```

## 2. Run a case

The cases under `test/` are complete: mesh, initial condition, boundary table,
`input.ini` and a `verify.py` that checks the result. The crossing-jets case with the
isotropic Gaussian closure is a good first one — it takes about a minute on one thread.

```bash
cd test/Doisneau/IG
./ICE.sh solve
```

`ICE.sh` creates `OUTPUT/`, copies the current `bin/ICE` into the case, and runs it.
`-p N` gives it `N` OpenMP threads and `-m N` runs it on `N` MPI ranks:

```bash
./ICE.sh -p 4 solve
```

Running `bin/ICE` directly from the case directory works too, as long as `OUTPUT/`
exists.

## 3. Watch it run

```
 Local time step 0-way coupled simulation

 ICE phase model:
 - IG particle families  -->    1

 ICE numerical scheme:
 - Space   --> MUSCL-SD with VANLEER flux limiter
 - Time    --> Explicit Runge-Kutta 2

 Boundary Conditions:
   Symmetry                       80
   Inflow                         20
   Extrapolation                  300

ICE  | Iter =       10 | Global iter =       10 | Density residual = 0.162364E-01
```

The first line says how the step is taken (local, so steady state) and how the phases
are coupled — `0-way` because this case has no gas field. The progress line then
reports the density residual in steady-state mode, and the simulation time and step in
time-accurate mode; `shell-diter` sets how often it appears.

## 4. Check the result

```bash
python3 verify.py
```

The case compares its output against a stored reference and prints `PASS` or `FAIL`.
The whole suite runs from the build directory:

```bash
ctest --test-dir build -j 5 --output-on-failure
```

## 5. Look at the output

`OUTPUT/part-field.tec` holds the solution: node coordinates, then the primitive
variables of every family, cell-centred. Open it in Tecplot, or switch `sol-format` to
`vtk binary` in `input.ini` for ParaView and VisIt.

`OUTPUT/part-residual-history.dat` has one row per `res-diter` iterations, with the
iteration number, the time and the five residual norms.

## Next steps

* **[Using ICE](../user/using.md)** — case layout, running in parallel, restarting.
* **[Input File](../user/input.md)** — the sections of `input.ini`.
* **[Input Parameters](../user/registry.md)** — every parameter, default and rule.
* **[Verification & Validation](../vv/index.md)** — what the shipped cases demonstrate.
