---
title: Verification & Validation
---

# Verification & Validation

Test cases used to verify the numerical implementation of ICE, by comparison against
exact solutions and against ICE's own single-block results.

The [code-verification](verification.md) cases and the [Berthon Riemann
problems](berthon.md) compare ICE against solutions that come from outside ICE, so they
can fail on their first run. The other cases compare against exact free streaming or
against ICE's own single-block solution.

## Test suite

| Test | Dim | Models | Physics | Verification | Regression case |
|---|---|---|---|---|---|
| [Crossing Jets](crossing-jets.md) | 2D | MK, IG, AG | Two free-streaming particle jets crossing at 90° | Exact free-streaming solution | `Doisneau/MK`, `Doisneau/IG`, `Doisneau/AG` |
| [Chimera Overset](chimera.md) | 2D | IG | Crossing jets on two overlapping, non-matching blocks | Single-block solution; matched-resolution control | `Doisneau/IG-chimera` |
| [Multi-block and MPI](multiblock-mpi.md) | 2D | IG | Crossing jets on four blocks joined by connections | Single-block solution; identical results on 1 to 4 MPI ranks | `Doisneau/IG-split4` |
| [Berthon Riemann Problems](berthon.md) | 1D | AG | Shock, contact and rarefaction waves of the Gaussian closure | Analytical wave patterns | `Berthon/SCS`, `Berthon/RCS`, `Berthon/RCR` |
| [Code Verification](verification.md) | 1D | MK | Drag and heat relaxation, periodic advection, every drag and Nusselt law | Closed-form solutions and independent RK4 integration | `verification/A` … `verification/E` |

## Running the tests

Each case lives under `test/` with a `verify.py` script: the crossing-jets cases compare
the density field against the case's stored reference solution (L2 norm, tolerance
$10^{-4}$), the Berthon cases compare six fields against the analytical wave pattern
(L1 norm), and the code-verification cases build their own inputs from scratch.

```bash
ctest --test-dir build -j 5 --output-on-failure   # run and verify every case
ctest --test-dir build -R DoisneauIG -V           # a single case
```

See [Testing](../development/testing.md) for the tiers, the labels and the pre-push hook.

The figures on these pages are produced from the same runs:

```bash
cd test/Doisneau
python3 plot_vv.py                  # writes docs/vv/images/*.svg and prints the quoted metrics
```
