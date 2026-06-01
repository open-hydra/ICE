# Boundary Conditions

Boundary conditions for ICE are specified in the `bc/` directory. Each block face is assigned a boundary type in `input.ini`.

## Available types

| Type | Description |
|------|-------------|
| `connection` | Block-to-block interface — continuity of flux across the boundary |
| `null` | Inactive (degenerate) face — used for 2-D cases on faces 5/6 |
| `wall` | Impermeable wall — zero normal particle flux |
| `inlet` | Prescribed inflow state for volume fraction, velocity, and temperature |
| `outlet` | Extrapolated outflow — zero-gradient condition |
| `symmetry` | Symmetry plane — normal velocity component reflected |

## Specification format

Boundary conditions are declared in the `[bc]` section of `input.ini` (or in a separate `bc/bc.ini` file, depending on the case layout). Each block face (`face1` through `face6`) is assigned a type and, where required, additional parameters.

---
