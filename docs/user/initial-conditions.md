# Initial Conditions

ICE reads its initial state from a single file, `INPUT/part-ic.tec`, which carries the
**mesh and the solution together**: node coordinates as the first three variables, then
the primitive variables of every family as cell-centred data. There is no separate mesh
input, and no per-block or per-family file.

The reader is chosen by `ic-format`, so the file is Tecplot ASCII by default and
`INPUT/part-ic.vtm` (plus `INPUT/vtk/*.vts`) when `ic-format` names VTK.

## Variable order

The variables follow the families in order, and within a family they are the primitive
variables of its closure:

| Closure | Variables, in order |
|---|---|
| MK | `rho_p` `u_p` `v_p` `w_p` `T_p` `n_p` |
| IG | `rho_p` `u_p` `v_p` `w_p` `P_p` `T_p` `n_p` |
| AG | `rho_p` `u_p` `v_p` `w_p` `P11_p` `P12_p` `P13_p` `P22_p` `P23_p` `P33_p` `T_p` `n_p` |

For two families, family 1's block is followed by family 2's. ICE names them with the
family index appended — `rho_p1`, `u_p1`, … `rho_p2` — when it writes, but the reader
goes by **position, not by name**: a file whose variables are in the wrong order is
read without complaint and gives a wrong answer.

| Symbol | Meaning | Units |
|---|---|---|
| $\rho_p$ | bulk density: condensed mass per unit volume of mixture | kg/m³ |
| $u_p, v_p, w_p$ | particle velocity | m/s |
| $P$, $P_{ij}$ | velocity dispersion (IG: scalar, AG: symmetric tensor) | Pa |
| $T_p$ | particle temperature | K |
| $n_p$ | number density | m⁻³ |

$\rho_p$ and $n_p$ together fix the particle radius,
$R_p = \big(\tfrac{3}{4\pi}\rho_p/(n\rho_{al})\big)^{1/3}$; there is no radius variable.
A cloud of radius $R_p$ at bulk density $\rho_p$ therefore needs

$$
n = \frac{\rho_p}{\rho_{al}\,\tfrac43\pi R_p^3}.
$$

For the Gaussian closures a strictly positive dispersion is needed even where the cloud
is nominally monokinetic; the shipped cases use $10^{-6}$ as a floor, which is also
what ICE writes into an inlet ghost cell.

## File layout

A Tecplot ASCII file, one `ZONE` per block, `DATAPACKING=BLOCK`, with the coordinates
nodal and the rest cell-centred:

```
VARIABLES ="x" "y" "z" "rho_p1" "u_p1" "v_p1" "w_p1" "T_p1" "n_p1"
ZONE T = Block1, I=101, J=101, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED)
<I*J*K x values>
<I*J*K y values>
<I*J*K z values>
<(I-1)*(J-1)*(K-1) rho_p1 values>
...
```

`I`, `J`, `K` count **nodes**, so a block of $n_x \times n_y \times n_z$ cells has
`I = nx+1` and so on. A 2-D case is written with `K=2`, a 1-D case with `J=K=2`; ICE
infers the dimensionality from the resulting cell counts.

!!! warning "The `VARIABLES` line must be one line, with the names separated"
    Writing it with Fortran list-directed output wraps it at 80 columns and runs the
    quoted names together. The reader then sees a single very long variable name, finds
    no blocks, and ICE fails afterwards in a way that gives no hint of the cause. Write
    it with an explicit `'(A)'` format and put a space between the names.

## Restarting

With `newrun = false` the initial state is read from `OUTPUT/part-field<ext>` instead,
and the simulation time is taken from the zone header. The format and variable order
are identical — it is the same writer.

## Property table

`INPUT/part-properties.dat` is optional and gives $\rho_{al}(T)$ and $c_s(T)$ for the
condensed material. It is a Tecplot point file with four columns — temperature, $c_p$,
density, enthalpy — indexed by **integer** temperature: row $i$ must be $T = i$, and
the table covers $1 \le T \le N$, saturating outside. If it is absent, the constants
`density` and `specific-heat` from `[ICE-Physics]` are used; ICE prints which of the two applies.

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"
ZONE T="A"
I=5000, F=POINT
   1.0   900.0  2700.0   ...
   2.0   900.1  2700.0   ...
```

## Gas field

`INPUT/gas.tec` is read only if it exists, and its presence is what turns on the source
terms. It must match the mesh block for block and supply ten variables per cell:
density, the three velocity components, temperature, $R$, $\gamma$, thermal
conductivity, viscosity and a time step. ICE reads it once and never modifies it.
