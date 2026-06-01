# Input File

The main configuration file for ICE is `input.ini`, an INI-format file parsed by [FiNeR](https://github.com/szaghi/FiNeR).

## File structure

`input.ini` is organised into named sections, each controlling a specific aspect of the solver.

```ini
[simulation]
  ...

[io]
  ...

[numerics]
  ...

[part-1]
  ...

[part-2]
  ...
```

## Section reference

For the full auto-generated parameter reference see [Input Parameters](registry.md).

### `[simulation]`

Global simulation parameters: time limits, restart options, parallelism settings.

### `[io]`

Output format and frequency: solution files (`sol_format`, `sol_diter`), backup files (`bck_format`, `bck_diter`), and probe output.

### `[numerics]`

Solver settings: Runge–Kutta stages, CFL number, time-step control mode.

### `[part-N]`

One section per particle group (named `part-1`, `part-2`, …). Specifies material properties (`rho_al`, `cs_al`, `lv_al`, `q_al`), heat-transfer model, and drag model for the group.

---
