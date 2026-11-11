---
title: Verification & Validation
---

# Verification & Validation

These pages ask two different questions, and it helps to keep them apart:

- **Verification**: does ICE solve its equations correctly? The case compares against a
  solution that does not come from ICE (a closed form, an analytical wave pattern, or an
  integration done independently in Python). It can fail on its first run, and a failure
  means something is wrong.
- **Regression**: does ICE still give the answer it gave before? The case compares against a
  stored solution of its own. It cannot tell that the answer is right, only that it has not
  changed.

Most cases do both: they establish the answer against an exact solution once, and keep a
stored reference so that a change is caught.

## The pages

The pages fall into three groups.

### Code verification

The models and the schemes, one at a time, against closed forms and independent integrations.
Every case builds its own inputs; see [Code Verification](verification/index.md) for how they
work and for the full list.

| Page | Cases | What it checks |
|---|---|---|
| [Drag and Heat](verification/exchange.md) | A, B, D, E, J | The drag and heat sources, every drag and Nusselt correlation, the time scheme's order, the step at the relaxation time |
| [Evaporation and Boiling](verification/phase-change.md) | F, M | Every evaporation model, the latent sink, the boiling clamp, the temperature inside the drop, the boiling branch |
| [Solidification](verification/solidification.md) | N, O | Supercooling, nucleation, the plateau, melting; a freezing front in a steady stream |
| [Combustion](verification/combustion.md) | P | Drops burning on the $d^2$ law down to their residue, ignition, the heat to the gas |
| [Material Properties](verification/materials.md) | K, L | The property table, and several materials side by side |
| [Transport](verification/transport.md) | C, G, H, I | Transport alone and with drag, a straining gas, a vortex; the order of the space scheme |

### Closures and waves

How the MK, IG and AG closures behave where they differ: crossing streams, shocks, and the
pressure tensor.

| Page | Dim | Closures | What it checks | Compared against | ctest |
|---|---|---|---|---|---|
| [Crossing Jets](crossing-jets.md) | 2D | MK, IG, AG | Free streaming, inlets, how each closure handles crossing streams | Exact free streaming (IG, AG), sticky merger (MK), stored reference | `DoisneauMK`, `DoisneauIG`, `DoisneauAG` |
| [Berthon Riemann Problems](berthon.md) | 1D | AG | Shocks, contacts and rarefactions of the Gaussian system | Analytical wave patterns | `BerthonSCS`, `BerthonRCS`, `BerthonRCR`, `BerthonHLLE` |
| [Gaussian Closures](gaussian-closures.md) | 1D, 2D | IG, AG | The pressure work in the energy flux, the directional wave speed, the reflected pressure tensor | A uniform temperature, the exact γ = 3 Riemann solution, momentum conservation | `ThermalIG`, `ThermalAG`, `RiemannAG`, `ReflectAG` |

### Grids and boundaries

The ways a block meets its neighbours and the rest of the domain, and the parallel
decomposition.

| Page | Dim | Closures | What it checks | Compared against | ctest |
|---|---|---|---|---|---|
| [Axisymmetric Wedge](axisymmetry.md) | 2D-axi | MK, IG, AG | The wedge side faces and the axis face, and the hoop pressure | Exact stationary states | `AxisMK`, `AxisIG`, `AxisAG` |
| [Chimera Overset](chimera.md) | 2D | IG | Overset interpolation between overlapping, non-matching blocks | The same case on one block, and a matched-resolution control | `DoisneauChimera` |
| [Multi-block and MPI](multiblock-mpi.md) | 2D | IG | Block connections and the rank decomposition | The single block; 1 to 4 ranks against each other | `DoisneauSplit4` |

## What is covered

Read together, the pages cover:

- every closure;
- all three source terms (drag, heat, mass) and every one of their correlations: twelve drag
  laws, seven Nusselt laws and five evaporation models;
- the drop's internal temperature, boiling, solidification and direct combustion;
- the transport operator and its order of accuracy, and the time integrator and its order;
- every way a block can talk to another one.

### What is not covered

- **Restarts, probes, grid sequencing, 3-D meshes and several families on a mesh** are checked
  by the implementation tests ([Testing](../development/testing.md)): a restarted run must
  reproduce the uninterrupted one byte for byte, a 3-D run on several threads or ranks the
  serial one, and so on. Those tests show that the code is consistent, not that its answer is
  right; no case on these pages measures their accuracy against an outside answer.
- **Implicit residual smoothing** has no case at all.
- **Evaporation** is verified only at zero slip, where the Sherwood and Nusselt corrections are
  inactive: the convective branch of each model is not verified.
- **Combustion in a coupled gas**: ICE fills the `heat` field a two-way coupled gas adds to its
  energy, and the unit test checks it, but no case here runs a gas that receives it.

## Running them

```bash
ctest --test-dir build -j 5 --output-on-failure   # everything, a few minutes
ctest --test-dir build -L verification            # only the cases with an outside answer
ctest --test-dir build -R DoisneauIG -V           # one case, with live output
```

The cases under `test/Doisneau/`, `test/Berthon/` and the other case folders run in their own
directory through `ICE.sh` and check themselves with a `verify.py`. The cases under
`test/verification/` generate their inputs from scratch into a scratch directory. See
[Testing](../development/testing.md) for the tiers, the labels and the pre-push hook.

## Regenerating the figures

```bash
cd test/Doisneau     && python3 plot_vv.py   # crossing jets and chimera, from the case outputs
cd test/Berthon      && python3 plot_vv.py   # the Riemann problems, from the case outputs
cd test/verification && python3 plot_vv.py   # code verification: runs what it needs, a few minutes
```

All three write into `docs/vv/images/` and print the numbers quoted in the text. The first two
read the outputs of the cases, so run the cases first.
