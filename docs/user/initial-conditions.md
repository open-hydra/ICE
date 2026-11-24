# Initial Conditions

ICE reads its initial state from one file, `INPUT/part-ic.tec`, which holds the
**mesh and the state together**: the node coordinates first, then the variables of
every family at the cell centres. There is no separate mesh file.

The format is set by `ic-format`: Tecplot ASCII by default, binary Tecplot
(`part-ic.szplt`), or VTK (`part-ic.vtm` with its block files).

## Variable order

The families follow one another, and each brings the variables of its closure:

| Closure | Variables, in order |
|---|---|
| MK | `rho_p` `u_p` `v_p` `w_p` `T_p` `n_p` |
| IG | `rho_p` `u_p` `v_p` `w_p` `P_p` `T_p` `n_p` |
| AG | `rho_p` `u_p` `v_p` `w_p` `P11_p` `P12_p` `P13_p` `P22_p` `P23_p` `P33_p` `T_p` `n_p` |

A family of a [solidifying](../theory/solidification.md#solidification) material has two more
variables after `n_p`: `f_p`, the frozen fraction, and `chi_p`, the nucleated fraction.
They may be left out of the file for every solidifying family at once: ICE then takes
a particle at or below `T-nuc` as solid, and above it as liquid.

!!! warning "Order, not names"
    ICE writes the names with the family number appended (`rho_p1`, `u_p1`, …, then
    `rho_p2`), but it **reads by position**. A file with the variables in another
    order is read without complaint, and gives a wrong answer.

| Symbol | Meaning | Units |
|---|---|---|
| $\rho_p$ | bulk density: particle mass per unit volume of mixture | kg/m³ |
| $u_p, v_p, w_p$ | particle velocity | m/s |
| $P$, $P_{ij}$ | velocity dispersion (IG: a scalar, AG: a symmetric tensor) | Pa |
| $T_p$ | particle temperature | K |
| $n_p$ | number density | m⁻³ |

There is no radius variable: $\rho_p$ and $n_p$ fix it. For particles of radius $R_p$
at bulk density $\rho_p$, set

$$
n_p = \frac{\rho_p}{\rho_{al}\,\tfrac43\pi R_p^3},
$$

with $\rho_{al}$ the density of the material.

The Gaussian closures need a strictly positive dispersion, even for a cloud meant to
be monokinetic. The shipped cases use $10^{-6}$, the value ICE also gives the
particles an inlet injects.

## File layout

A Tecplot ASCII file has one `ZONE` per block, in `BLOCK` packing, with the
coordinates at the nodes and the rest at the cell centres:

```
VARIABLES ="x" "y" "z" "rho_p1" "u_p1" "v_p1" "w_p1" "T_p1" "n_p1"
ZONE T = Block1, I=101, J=101, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED)
<I*J*K x values>
<I*J*K y values>
<I*J*K z values>
<(I-1)*(J-1)*(K-1) rho_p1 values>
...
```

`I`, `J` and `K` count **nodes**: a block of $n_x \times n_y \times n_z$ cells has
`I = nx+1`, and so on. A 2-D case has `K=2`, a 1-D case `J=K=2`; ICE infers the
dimensions from the cell counts.

!!! warning "Write the `VARIABLES` line as one line, with spaces between the names"
    Fortran list-directed output wraps it at 80 columns and runs the quoted names
    together. ICE then sees one very long name, finds no blocks, and fails later with
    no hint of why. Write it with an explicit `'(A)'` format.

A restart reads the same layout from `OUTPUT/part-field.*`: see [Restarting](using.md#restarting).

## Gas field

`INPUT/gas.tec` is read only if it exists, and its presence switches the exchange
with the gas on. It must match the mesh block by block, with ten variables per cell:
density, the three velocity components, temperature, $R$, $\gamma$, thermal
conductivity, viscosity and a time step. ICE reads it once and never changes it.
