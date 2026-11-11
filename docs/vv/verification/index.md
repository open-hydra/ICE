# Code Verification

These cases compare ICE against answers that do not come from ICE: a closed form, or an
accurate integration of the equations ICE is meant to solve. They can fail on their first run,
and a failure means that a model or a scheme is wrong.

## How the cases work

- Each case is a `run.py` under `test/verification/<letter>-<name>/`.
- It writes its own mesh, initial condition, boundary conditions and `input.ini` into a
  scratch directory, runs ICE there and compares. The repository holds no case data and no
  reference to regenerate.
- Most cases use the MK closure; the solidification cases also run IG and AG.
- Each case states its tolerance and why it was chosen. Where the pages list *what fails it*,
  that is the error a broken model would make, which is how the tolerance was set.

The problems are of two kinds:

- **A uniform cloud** (A, B, D to F, J, K, M, N): nothing varies in space, so the exact answer
  is the same in every cell and the mesh only exercises the flux loops. These cases measure
  the source terms.
- **A cloud carried by a flow** (C, G to I, O): a localized cloud crosses a frozen gas field,
  or a stream reaches a steady state, so the answer varies from cell to cell and the mesh
  resolution is part of what is measured.

L runs two or three families side by side, each on its own material, some at rest and some
fed by an inlet.

## The cases

Sorted by letter; each page groups its cases by what they check.

| Case | ctest | What it pins down | Page |
|---|---|---|---|
| A | `DragRelaxation` | The momentum source, and the order of the time scheme | [Drag and Heat](exchange.md#a-relaxation-under-stokes-drag) |
| B | `ThermalRelaxation` | The energy source at $Nu = 2$ | [Drag and Heat](exchange.md#b-thermal-relaxation) |
| C | `SineAdvection` | Transport alone, and the order of the space scheme | [Transport](transport.md#c-sinusoidal-advection-on-a-periodic-mesh) |
| D | `DragMatrix` | All twelve drag laws, up to $Ma = 2.6$ | [Drag and Heat](exchange.md#d-every-drag-correlation) |
| E | `HeatMatrix` | All seven heat laws | [Drag and Heat](exchange.md#e-every-nusselt-correlation) |
| F | `Evaporation` | All five evaporation models, the latent sink, the boiling clamp, the table's `Psat`, `Lv`, `Mv`, `Cpv` | [Evaporation and Boiling](phase-change.md#f-evaporation) |
| G | `CloudTranslation` | Transport and drag together | [Transport](transport.md#g-a-cloud-released-into-a-uniform-gas) |
| H | `LinearStrain` | A non-trivial particle velocity field, and the small-Stokes limit | [Transport](transport.md#h-a-cloud-in-a-straining-gas) |
| I | `VortexCloud` | Two dimensions, over four decades of Stokes number | [Transport](transport.md#i-a-cloud-in-a-prescribed-vortex) |
| J | `TauLimit` | `tau-factor`: the step bounded by $\tau_p$; nearly empty cells left out | [Drag and Heat](exchange.md#j-the-step-at-the-relaxation-time) |
| K | `PropertiesTable` | The property table: interpolation, its ends, a varying specific heat | [Material Properties](materials.md#k-the-property-table) |
| L | `Materials` | Each family on its own material: table zone, model tokens, inlet records | [Material Properties](materials.md#l-several-materials) |
| M | `DropHeating` | The temperature inside the drop (`P2T`), the boiling branch (`ZGR`) | [Evaporation and Boiling](phase-change.md#m-drop-heating-and-boiling) |
| N | `SolidCooling` | Supercooling, nucleation, the plateau, the solid, melting | [Solidification](solidification.md#n-a-closed-cell) |
| O | `SolidAdvection` | The nucleation front and the plateau in a steady stream | [Solidification](solidification.md#o-a-steady-stream) |

The letters are the directory names under `test/verification/`, given in the order the cases
were written.

## Running them

```bash
ctest --test-dir build -L verification --output-on-failure   # these and the other cases with an outside answer
ctest --test-dir build -R DropHeating -V                     # one case, with live output
```

- A, B, C, J, K, L, M, N and O also carry the `fast` label, so the pre-push hook runs them.
- The rest are in the default tier only: D, E and F because each correlation is a separate
  run of the solver, and G, H and I because each refinement level is. I, in two dimensions, is
  the longest.
- `ICE_KEEP_WORK=1` keeps the generated case directories instead of deleting them on
  success: that is how to look at what a failing case actually produced.

## The figures

```bash
cd test/verification && python3 plot_vv.py [figure ...]
```

| Figure | Page |
|---|---|
| `relaxation`, `correlations` | [Drag and Heat](exchange.md) |
| `evaporation` | [Evaporation and Boiling](phase-change.md) |
| `translation`, `strain`, `vortex`, `orders` | [Transport](transport.md) |

With no argument the script draws all of them, in a few minutes, most of it in the vortex. It
does not read the output of a ctest run: these cases build their inputs on the fly and remove
them again, so there is nothing left to plot. It runs what each figure needs itself, and
prints the numbers quoted on the pages.
