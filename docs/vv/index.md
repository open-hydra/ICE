---
title: Verification & Validation
---

# Verification & Validation

Test cases used to verify the numerical implementation of ICE, by comparison against
exact solutions and against ICE's own single-block results.

## Test suite

| Test | Dim | Models | Physics | Verification | Regression case |
|---|---|---|---|---|---|
| [Crossing Jets](crossing-jets.md) | 2D | MK, IG, AG | Two free-streaming particle jets crossing at 90° | Exact free-streaming solution | `Doisneau/MK`, `Doisneau/IG`, `Doisneau/AG` |
| [Chimera Overset](chimera.md) | 2D | IG | Crossing jets on two overlapping, non-matching blocks | Single-block solution; matched-resolution control | `Doisneau/IG-chimera` |

## Running the tests

Each case lives under `test/` with a `verify.py` script that compares the density field
against the case's stored reference solution (L2 norm, tolerance $10^{-4}$).

```bash
cd test
./test.sh -p 8 check all            # run and verify every case with 8 OpenMP threads
./test.sh -p 8 check Doisneau/IG    # a single case
```

The figures on these pages are produced from the same runs:

```bash
cd test/Doisneau
python3 plot_vv.py                  # writes docs/vv/images/*.svg and prints the quoted metrics
```
