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

A family of a [solidifying](../theory/physics.md#solidification) material has two more
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

## Property table

`INPUT/part-properties.dat` gives the properties of the materials as functions of
temperature. It is optional with one material, and required with several. ATLAS writes
it; this is its layout:

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"
ZONE T="A"
I=5000, F=POINT
   1.0   900.0  2700.0    900.0
   2.0   900.0  2700.0   1800.0
   ...
```

- **One zone per material**, in the order of the [phase file](input.md#materials),
  all on the same temperatures.
- **The rows sit on consecutive integer kelvins**, at least two of them, from any
  first temperature. Between the rows the values are linear in $T$; outside the table
  they keep the end values.
- **Packing:** `F=POINT` (one row per line, as above) or `BLOCK` (one column after
  the other). A zone that gives neither is `BLOCK`, as Tecplot reads it.

### Columns

`Temperature` comes first; the others may be in any order.

| Column | Unit | Required | Used for |
|---|---|---|---|
| `Temperature` | K | yes | the rows |
| `Cp` | J/(kg K) | yes | the specific heat of the material |
| `Density` | kg/m³ | yes | the density of the material |
| `Enthalpy` or `Enthalpy_abs` | J/kg | one of the two | the energy of the particles (below) |
| `Psat` | Pa | no | the saturation pressure, in place of the Clausius-Clapeyron curve |
| `Lv` | J/kg | no | the latent heat, against temperature |
| `Mv` | kg/kmol | no | the vapour molar mass, the same on every row |
| `Cpv` | J/(kg K) | no | the vapour specific heat, against temperature |
| `Emissivity` | – | no | the emissivity, against temperature |
| `Conductivity`, `Viscosity` | W/(m K), Pa s | no | the liquid conductivity and viscosity of a `P2T` material |
| `Sigma` | N/m | no | recognised and not used |

Other columns are ignored. `Enthalpy` is relative ($c_p T$ for a constant $c_p$);
`Enthalpy_abs` is absolute, and its offset is kept as the material's datum. When $c_p$
varies, the particle energy is the enthalpy column itself, so its datum does not matter
([Energy and temperature](../theory/governing-equations.md#energy-and-temperature)).

The optional columns are read only where they apply: `Psat`, `Lv`, `Mv` and `Cpv` for
a material that evaporates, `Conductivity` and `Viscosity` for one that runs
`liquid-conduction = P2T`, `Emissivity` for every material. A column of zeros counts
as absent. Emissivities above 1 are cut to 1, with a warning.

### Columns and INI keys together

- **`density` and `specific-heat`** may be left out with a table. If given, a value
  must equal the material's constant column; against a column that varies it stops the
  run. Without a table, the one material takes these two constants. ICE prints which
  of the two it uses.
- **`latent-heat`, `vapour-molar-mass`, `vapour-specific-heat` and `emissivity`** must
  be 0 for a material whose table has the same quantity as a column, even if the two
  agree. With two materials, `latent-heat = 0 2.26e6` lets the table set A's and gives
  B, which has no column, its value.
- **`k-liq` and `mu-liq`** on the material line must be absent (or 0) when the table
  has the `Conductivity` and `Viscosity` columns.

At setup ICE prints the tabulated vapour values of each evaporating material: $L_v$ at
298.15 K, which becomes its scalar latent heat, and $L_v$, $M_v$ and $c_{p,v}$ at the
boiling temperature.

### Checks

ICE refuses the table, with the reason, when:

- **ORION cannot read the file**: a row with a number too many or too few, a zone that
  holds another number of rows than it announces, anything after the last row. ORION's
  message says what and where;
- **the columns** are missing one, name one twice, or name both enthalpies;
- **the zones** are not one per material, are not one line of rows (`J=` or `K=` above
  1), or are not on the same temperatures;
- **the rows** go below 0 K or are off the integer kelvins (by more than $10^{-6}$ K);
- **the values** are not physical: a density or $c_p$ that is not positive, an
  enthalpy that does not increase with $T$, or disagrees with $c_p$ (with a constant
  $c_p$, $h = c_p T + h_{off}$ on every row; with a varying one, the trapezoidal
  integral of $c_p$ within 0.1 % per step), or a relative `Enthalpy` with an offset;
- **`Psat`**, for an evaporating material, is not finite, is negative, constant or
  decreasing, or does not give `boiling-pressure` within a factor 2 at
  `boiling-temperature` (which must lie inside the table);
- **`Lv`, `Mv`, `Cpv`** are not finite, negative, zero on some rows only, given without
  `Psat`, an `Mv` that varies, or an `Lv` whose Clausius-Clapeyron slope at the boiling
  temperature is not within a factor 2 of the `Psat` column's. An `Lv` column without
  `Mv` needs `vapour-molar-mass`, which the slope uses;
- **`Emissivity`, `Conductivity`, `Viscosity`** are zero on some rows only, or not
  positive.
