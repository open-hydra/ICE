# Material Properties

Two cases on where ICE takes a material's properties from. K reads them from the property
table, and checks the interpolation and the energy of a varying specific heat. L runs several
materials side by side, so that a family reading another family's material, model or inlet
record shows directly in its own field.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [K](#k-the-property-table) | `PropertiesTable` | Linear interpolation between the table's rows, a table that starts above 1 K, saturation past its ends, the energy of a varying specific heat | Closed form and RK4 |
| [L](#l-several-materials) | `Materials` | Each family on its own material: its table zone, its model tokens, its inlet records, its inlet density; the `401` inlet of several families | Closed form and case F's |

## K. The property table

### At a fixed temperature

**Setup.** Case A's cloud, with its density taken from `INPUT/part-properties.dat` instead of
the INI. With `heat-transfer = NoHeat` and no radiation nothing heats the particles, and the
drag work cancels the kinetic energy it produces. So $T_p$ keeps its initial value, and the
table is read at one temperature for the whole run.

**What is measured.** The velocity is case A's exponential with
$\tau_p = \rho_{al}(T_p)\,d_p^2/(18\mu_g)$, and $\rho_{al}(T_p)$ is what the case measures.
ICE recovers the radius from $\rho_p$ and $n$ through the table's density, so at fixed
$(\rho_p, n)$ the relaxation time goes as $\rho_{al}^{1/3}$.

| Leg | Table | $T_p$ | $\rho_{al}(T_p)$ | Relative error on $u_p(\tau_p)$ |
|---|---|---|---|---|
| K1 | 1 to 5000 K; $\rho$ = 2000 up to 300 K, 1000 from 301 K | 300.4 K | 1600 (linear) | $9.8\times10^{-8}$ |
| K1c | the same rows, $\rho$ = 1600 everywhere | 300.4 K | 1600 | $9.7\times10^{-8}$ |
| K6 | 280 to 400 K, constant | 300.4 K | 1500 | $9.7\times10^{-8}$ |
| K7 | 280 to 400 K, $\rho = 1500 + 5\,(T - 280)$ | 250 K | 1500 (the end value) | $9.7\times10^{-8}$ |

- **The tolerance** is $10^{-4}$. RK2 at $\Delta t = \tau_p/1000$ leaves
  $e^{-1}(\Delta t/\tau_p)^2/6 = 6.1\times10^{-8}\,u_g$ on $u_p(\tau_p)$, that is
  $9.7\times10^{-8}$ of $u_p(\tau_p)$ itself, which is what the four legs measure.
- **What fails it.** The nearest row would give K1 $\rho_{al} = 2000$ and put $u_p(\tau_p)$
  4.3 % low. A linear extrapolation below $T_{min}$ would put K7 2.0 % high, and the row at
  $T_{max}$ 6.5 % low.
- **K1c** is K1's control: the same run with nothing to interpolate.
- **The temperature** must stay within $10^{-4}$ K of its initial value, well under the 0.1 K
  that would change the nearest row. It moves by round-off.

### A specific heat that varies

The energy of a material whose `Cp` column varies is $\rho_p e(T_p)$ with $e = h - h_{off}$
([Energy and temperature](../../theory/governing-equations.md#energy-and-temperature)). The
tables below have $c_s = 1000 + 2T$ on the rows 1 to 5000 K, $h$ its exact integral (whose
values are integers, so that every subtraction is exact), and a density of 1000 kg/m³.

| Leg | Setup | Check | Measured |
|---|---|---|---|
| K2 | No gas; a uniform cloud at 10 m/s and 500.37 K, 200 RK2 steps | $\lvert T_p - T_{p0}\rvert \le 10^{-9}$ K | 0 |
| K2c | K2 with $c_s$ = 2000 on every row | the same | 0 |
| K3 | The cloud moving with a gas at 800 K, $T_{p0}$ = 300 K, $Nu = 2$, no drag, to $t = \tau_T$ | $T_p$ within $10^{-6}$ of the 500 K gap of the reference | $4.7\times10^{-8}$ |
| K3c | K3 with $c_s$ = 2000 on every row | the same | $6.1\times10^{-8}$ |
| K4 | K3's table as `Enthalpy_abs`, $h - 1.5\times10^7$ J/kg | `part-field.tec` bit for bit K3's | identical |
| K5 | K3's table with its columns in the order `Temperature`, `Enthalpy`, `Density`, `Cp` | the same | identical |

- **The reference** of K3 integrates $\rho_p\,de/dt = 4\pi k_g R_p\,n\,(T_g - T(e))$ with RK4
  in $e$, $T(e)$ the inverse of the table's line, and $\tau_T = \rho_{al} c_s d_p^2/(12 k_g)$
  with $c_s$ = 2000.
- **What fails it.** An energy $\rho_p c_s(T_p)\,T_p$ is what these legs detect: its inversion
  is not a round trip, so K2's cloud drifts 201 K towards 299 K in 200 steps, and K3's never
  heats (error 0.64 of the gap). K2c and K3c are the controls with nothing to vary.
- **K2's tolerance.** A round trip $T \to e \to T$ returns $T$ to a few ulp, $10^{-13}$ K at
  500 K, and the 400 of the run stay below $10^{-10}$ K. The output holds 17 digits, which give
  the double back, so the check resolves its last bit ($10^{-13}$ K).
- **K3's tolerance.** RK2 at $\Delta t = \tau_T/1000$ leaves
  $e^{-1}(\Delta t/\tau_T)^2/6 = 6\times10^{-8}$ of the gap at $t = \tau_T$, which the
  constant-$c_s$ control measures too. The kinks of $e$ at the rows add at most $10^{-9}$ each
  over the 300 rows crossed. $10^{-6}$ sits a decade above both.

## L. Several materials

The phase file names the materials, one line each, `<name> <groups> [key=value ...]`, and the
families map onto them in that order; the property table gives one zone per material. Every
leg runs two or three MK families side by side in one case.

The legs come in five groups:

| Group | What it checks |
|---|---|
| [L1](#properties) | Each family heats with its own material's properties |
| [L2](#model-tokens) | The per-material model tokens, and the INI defaults |
| [L3, L4](#inlet-records) | Each family reads its own inlet records and its own inlet density |
| [L5](#latent-heat) | Each material evaporates with its own latent heat |
| [L6](#the-401-inlet-of-several-families) | The `401` inlet of several families: the mixture rule |

### Properties

Families at rest heat in a gas at 400 K ($Nu = 2$). Each must follow
$T_g + (T_0 - T_g)\,e^{-t/\tau}$ with its own $\tau = \rho_p c_s/(4\pi k_g R_p n)$, and $R_p$
from its own density.

| Leg | Phase file | What it checks | Measured |
|---|---|---|---|
| L1 | `A 1` / `B 1`, zones $c_s$ = 1000, $\rho$ = 2000 and $c_s$ = 2000, $\rho$ = 1000 | Each family on its own $\tau$ | $6.1\times10^{-8}$, $2.2\times10^{-8}$ of the gap |
| L1b | `A 2` / `B 1`, three families | Families 1 and 2 heat as A, family 3 as B | the same |
| L1c | `A 2`, one zone | Control: two families on one material | the same |

The tolerance is $10^{-6}$ of the initial gap. RK2 at $\Delta t = \tau_A/1000$ leaves at most
$(t/\tau)(\Delta t/\tau)^2/6 = 6\times10^{-8}$ of it, and a family that heats with the other
material misses by 0.165.

### Model tokens

| Leg | Phase file | What it checks | Measured |
|---|---|---|---|
| L2 | `A 1 evaporation=CEM` / `B 1`, no `[ICE-Physics] evaporation` | A evaporates as case F's CEM (the $d^2$ slope); B keeps $\rho_p$ and $n$ bit for bit | slope $1.6\times10^{-8}$ |
| L2b | L2 with `evaporation = none` in the INI | The token overrides the INI default: `part-field.tec` bit for bit L2's | identical |
| L2d | `A 1` / `B 1 evaporation=none`, `evaporation = CEM` in the INI | The INI default applies to A, the token switches B off | slope $1.6\times10^{-8}$ |
| L2c, L2e | `A 1 evaporation=CEM`; `A 1` with `evaporation = CEM` | Controls: one material, the token alone and the INI alone | slope $1.6\times10^{-8}$ |
| L2f, L2g | `A 1 evaporation=CEM` / `B 1 evaporation=CEM`, two zones of one material, `boiling-temperature = 373.15 403.15` | One boiling temperature per material; L2g gives the same values under the alias `Tboil`: `part-field.tec` bit for bit L2f's | identical |
| L2h | L2f with the two values swapped | Control: the field changes, so each value reaches its own material | differs |

The slopes take case F's tolerance, $10^{-6}$.

### Inlet records

| Leg | Phase file | What it checks | Measured |
|---|---|---|---|
| L3 | `A 2`, face-1 inlets in ATLAS's order (the records of family 1, radius 10 µm, then family 2, 20 µm) | $n_1/n_2 = (r_2/r_1)^3 = 8$ in the cells the inlet fills; each family's $n = \rho_p/(\rho_{mat}\tfrac43\pi r^3)$ | ratio exact; $n$ to $8.9\times10^{-16}$ |
| L3b | L3 on two mesh blocks (mesh block, then family, then faces), radii 10/20 and 20/10 µm | 8 in block 1, 1/8 in block 2 | ratio exact and $2.2\times10^{-16}$; $n$ to $1.1\times10^{-15}$ |
| L3c, L3d | Controls: equal radii; one record per face on two blocks | Ratio 1 | exact |
| L4 | `A 1` / `B 1`, zones $\rho$ = 2000 and 1000, one record per face | Each family's inlet $n$ uses its own table density: $n_1/n_2 = 1/2$ | ratio exact; $n$ to $8.9\times10^{-16}$ |
| L4c | `A 2`, one zone | Control: ratio 1 | exact |

L3 and L4 compare ratios that are exact in binary ($r_2 = 2r_1$, $\rho_A = 2\rho_B$), through a
flux and a limiter that are homogeneous of degree one in $n$, so $10^{-12}$ bounds their
round-off. A cell counts as filled at half the inlet density.

### Latent heat

| Leg | Phase file | What it checks | Measured |
|---|---|---|---|
| L5 | `A 1 evaporation=d2-law` / `B 1 evaporation=d2-law`; zone A with `Lv` = $2.0\times10^6$, `Mv` and `Cpv` columns, zone B with the three at zero; `latent-heat = 0 2.26e6` | A evaporates with the latent heat of its columns, B with its INI value: each follows case F's $d^2$ slope with its own $L_v$ ($K_A/K_B = 1.116$) | slope $9.5\times10^{-10}$, $7.5\times10^{-10}$ |

The slopes take case F's tolerance, $10^{-6}$.

### The `401` inlet of several families { #the-401-inlet-of-several-families }

**Setup of L6.** `A 2`, with face-1 `401` records in ATLAS's order (family 1 `krho` = 0.2,
family 2 0.1, $kV = kT = 1$), in a uniform gas. Each family starts with $10^{-4}\,krho_i$
kg m⁻³. Field 1 is written as ATLAS writes it, $krho_i/(1-\sum_j krho_j)$, and $G = 0$.

| Leg | Setup | What it checks | Measured |
|---|---|---|---|
| L6 | as above | Family $i$ enters at $krho_i/(1-\sum_j krho_j)$ of the gas density, the mixture rule: 0.2/0.7 and 0.1/0.7 in cell 1, $\rho_{p1}/\rho_{p2} = 2$ in every cell | cell 1 $3.2\times10^{-9}$; ratio exact |
| L6c | `A 1`, `krho` = 0.2 | Control: one family, 0.25 in cell 1, where the two rules agree | $3.2\times10^{-9}$ |
| L6g | `A 1`, `krho` = 0.2 (field 1 0.25), $G = \rho_g U$ = 12 kg s⁻¹ m⁻² | The family takes field 1 times $(\rho_g U + G)$ through the face: 0.5 of the gas density in cell 1 (0.25 without the $G$ term) | $3.2\times10^{-9}$ |
| L6z | L6g in a gas at rest | No normal speed: $G$ left out, not divided by zero; every density finite and `part-field.tec` bit for bit that of $G = 0$ | identical |
| L6r | L6 with `krho` = 0.6 and 0.5, and with 0.2 and −0.1; L6's records without their 10th field $G$ | A negative field 1 (a sum of 1.1 makes it negative), or a `401` line without $G$ (a file of an ATLAS before the mixture rule over the dispersed phases), stops the setup with the `401` inlet message | refused |

- **The ratio** is exact for the same reason as L3's. The content keeps every cell above the
  empty-cell density, $10^{-6}$ kg m⁻³, below which the reconstruction is first order (a switch
  that is not scale-free), and starts in the inlet's ratio, so $10^{-12}$ bounds it too.
- **At the end time** cell 1 is still settling: its value is $3.2\times10^{-9}$ off; at four times the end time
  it is $1.9\times10^{-15}$ off.
- **Cell 1** is judged at $10^{-3}$ relative. A family loaded as if alone,
  $krho_i/(1-krho_i)$, misses there by 12.5 % and 22.2 %. L6g is judged on the same bound,
  which the $G$ term left out misses by 50 %.
