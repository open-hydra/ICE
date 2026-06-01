# Quick Start

This guide walks you through running your first ICE simulation after a successful [installation](installation.md).

## 1. Verify the build

After building, the `bin/ICE` executable should be present:

```bash
ls bin/ICE
```

## 2. Prepare a test case

!!! note
    Test cases are located under `test/`. Each case directory contains the input files required to run ICE.

Navigate to a test case:

```bash
cd test/mono-disperse
```

## 3. Run the simulation

```bash
./ICE.sh solve
```

Or invoke the executable directly:

```bash
../bin/ICE
```

## 4. Inspect the output

Solution files are written in the format specified in `input.ini` (Tecplot `.plt` or VTK `.vtu`). Open them with your preferred post-processing tool (Tecplot, ParaView, VisIt).

## Next steps

* **[User Guide](../user/using.md)** – learn how to configure and run your own cases.
* **[Input File Reference](../user/input.md)** – full description of all `input.ini` parameters.

---
