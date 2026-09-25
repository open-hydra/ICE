---
title: Verification & Validation
---

# Verification & Validation

Two different questions are asked on these pages, and it is worth keeping them apart.

**Verification** asks whether ICE solves its equations correctly, by comparing against a
solution that does not come from ICE: a closed form, an analytical wave pattern, or an
integration performed independently in Python. These cases can fail on their first run,
and a failure means something is wrong.

**Regression** asks only whether ICE still gives the answer it used to, by comparing
against a stored solution of its own. These cases cannot tell you the answer is right —
only that it has not changed.

Most pages here do both: they establish the answer against an exact solution once, and
then keep a stored reference so a change is caught.

## Test suite

| Test | Dim | Closures | What it exercises | Compared against | Cases |
|---|---|---|---|---|---|
| [Crossing Jets](crossing-jets.md) | 2D | MK, IG, AG | Free streaming, inlets, how each closure handles crossing streams | Exact free-streaming solution | `Doisneau/MK`, `Doisneau/IG`, `Doisneau/AG` |
| [Chimera Overset](chimera.md) | 2D | IG | Overset interpolation between overlapping non-matching blocks | The same case on a single block, plus a matched-resolution control | `Doisneau/IG-chimera` |
| [Multi-block and MPI](multiblock-mpi.md) | 2D | IG | Block connections and the rank decomposition | The single block; and 1 to 4 ranks against each other | `Doisneau/IG-split4` |
| [Berthon Riemann Problems](berthon.md) | 1D | AG | Shocks, contacts and rarefactions of the Gaussian system | Analytical wave patterns | `Berthon/SCS`, `Berthon/RCS`, `Berthon/RCR` |
| [Axisymmetric Wedge](axisymmetry.md) | 2D-axi | MK, IG, AG | The wedge side faces and the axis face of an axisymmetric mesh, and the hoop pressure | Exact stationary states | `Axis/MK`, `Axis/IG`, `Axis/AG` |
| [Gaussian Closures](gaussian-closures.md) | 1D, 2D | IG, AG | The pressure work in the energy flux, the directional wave speed of Rusanov and of the time step, the reflected pressure tensor at a symmetry plane | A uniform temperature, the exact γ = 3 Riemann solution, momentum conservation | `Thermal/IG`, `Thermal/AG`, `Riemann/AG`, `Reflect/AG` |
| [Code Verification](verification.md) | 1D, 2D | MK | Drag and heat relaxation, periodic advection, every drag and Nusselt correlation, every evaporation model, the property table, several materials side by side, and clouds transported through uniform, straining and rotating carrier fields | Closed forms and independent RK4 integrations | `verification/A` … `verification/L` |

Read together they cover: every closure, all three source terms and every one of
their correlations — twelve drag laws, seven Nusselt laws and five evaporation models —
the transport operator and its order of accuracy, the time integrator and its order,
and every way a block can talk to another one.

### What is not covered

Restart, probes, grid sequencing and implicit residual smoothing have no case. Nor does
any three-dimensional configuration — every case here is 1-D or 2-D — and several
families run side by side only in the uniform clouds of `verification/L`. Evaporation is covered only at zero slip, where the Sherwood and
Nusselt corrections are inactive; the convective branch of each model is not verified.

## Running them

```bash
ctest --test-dir build -j 5 --output-on-failure   # everything, about a minute
ctest --test-dir build -L verification            # only the cases with an outside answer
ctest --test-dir build -R DoisneauIG -V           # one case, with live output
```

Each case under `test/Doisneau/` and `test/Berthon/` runs in its own directory through
`ICE.sh` and checks itself with a `verify.py`; the cases under `test/verification/`
generate their inputs from scratch into a scratch directory. See
[Testing](../development/testing.md) for the tiers, the labels and the pre-push hook.

## Regenerating the figures

The figures on these pages come from the same runs, after the cases have been run:

```bash
cd test/Doisneau && python3 plot_vv.py   # crossing jets and chimera
cd test/Berthon  && python3 plot_vv.py   # the Riemann problems
```

Both write into `docs/vv/images/` and print the numbers quoted in the text.
