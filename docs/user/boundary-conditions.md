# Boundary Conditions

ICE reads its boundary conditions from `INPUT/<phase>-bc.txt` (`INPUT/bc.txt` for an
unnamed phase), written by [ATLAS BCB](https://github.com/open-hydra/ATLAS). Boundary
types are chosen per block face in the ATLAS input; ICE only consumes the resulting file.

## File format

One record per boundary face cell, in ATLAS's standard layout: a header line

```
block  i  j  k  face  code
```

followed by zero or more payload lines depending on `code`. Each record applies to every
particle group. The full specification is in the ATLAS BCB output reference.

## Supported codes

| Code | Type | Payload | ICE treatment |
|------|------|---------|---------------|
| `0` | Null | none | Ignored (degenerate faces, e.g. faces 5/6 in 2-D) |
| `101` | Connection | `b2 i2 j2 k2 f2 c1 c2 c3 c4` | Ghost cells copied from the neighbouring block |
| `201` | Periodic | as `101` | As `101` |
| `102` | Chimera (overset) | donor counts per ghost layer, then `b i j k weight` per donor | Ghost cells set to the weighted conservative blend of their donor cells |
| `200` | Axisymmetry | none | Ignored |
| `300` | Symmetry | none | Normal velocity reflected; except on face 3, it switches to extrapolation depending on the direction of particle motion (see `compute_ghost`) |
| `400` | Extrapolation | none | Zero-gradient outflow |
| `401` | Inlet, loading ratios | `krho kV alpha beta kT rp ...` | Mass loading, velocity and temperature as ratios of the local gas |
| `402` | Inlet, absolute | `gp v alpha beta Tp rp ...` | Mass flux, speed and temperature prescribed |
| `403` | Inlet, speed ratio | `gp kV alpha beta Tp rp ...` | As `402`, with speed a ratio of the local gas speed |

Any other code stops ICE at setup. For inlets, `alpha`/`beta` written as `normal` inject
along the face normal. If particles are leaving through an inlet face, ICE switches it to
extrapolation; with zero mass flux, it acts as symmetry.

## Chimera

For `102` faces, ATLAS BCB finds the donor cells overlapping each of the two ghost layers
and writes their volume weights, normalised to 1 per layer. ICE stops at setup if any ghost
layer has no donors. Donor blocks must cover the ghost cells completely: ATLAS normalises
the weights over the donors it finds, so partial coverage is not reported (check
`chimera.log`).
