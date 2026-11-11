# Combustion

One case on drops that burn with `combustion = direct`: P, a closed cell where the drops burn alone, on the $d^2$
law, down to their residue. A unit test checks the cell source against the paper's burn rate with slip.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [P](#p-a-closed-cell) | `DirectCombustion` | The $d^2$ law at the burning temperature, the residue, ignition from below, the refusals | Closed form |
| [Unit](#the-cell-source) | `CombustionUnit` | The burn rate with slip, the heat to the gas, the inert source | Gallier and Godfroy's eq. 4 |

## The material

Gallier and Godfroy's aluminium (J. Propuls. Power 25(2), 2009, Tables 1 and 2):

| Property | Value |
|---|---|
| Drop diameter $d_0$, residue $d_{res}$ | 120 µm, 60 µm |
| Density $\rho_m$, specific heat $c_p$ | 2000 kg/m³, 1177 J/(kg K) |
| Heat of reaction $Q$, latent heat $L_v$ | 9.53 MJ/kg, 10.8 MJ/kg |
| Burning temperature $T_{sat}$ | 2791 K |
| Gas | 3450 K, 8.71 MPa; $\mu$ = 9.1×10⁻⁵ Pa s, $Pr$ = 0.4, $\gamma$ = 1.16, $c_{p,g}$ = 2000 J/(kg K); one species |

The phase line is `A 1 combustion=direct q-comb=9.53e6 d-res=6e-5`. With $Nu = 2$ and
$B = (c_{p,g}(T_g - T_{sat}) + Q)/L_v$ = 1.004444, every regime is a closed form:

- burning, $T_p = T_{sat}$ and $d^2 = d_0^2 - K t$ with $K = 8 k_g \ln(1 + B)/(c_{p,g}\rho_m)$ = 6.328×10⁻⁷ m²/s,
  down to $d_{res}$ at $t_{res}$ = 17.07 ms;
- without a residue the drop would burn out at $d_0^2/K$ = 22.8 ms, the 23 ms of the paper's eq. 10;
- the residue relaxes on $\tau_{res} = m_{res} c_p/(\pi d_{res} k_g Nu)$ = 1.552 ms;
- a drop below $T_{sat}$ heats on $\tau_0$ = 6.208 ms and ignites when it reaches $T_{sat}$.

## P. A closed cell

**Setup.** Eight uniform cells moving with the gas (`drag = NoDrag`), so each integrates the source alone;
`heat-transfer = Stokes`, emissivity 0, RK2 with $\Delta t = 10^{-5}$ s. The gas phase file `INPUT/phase.txt`
declares one species. The closed form is taken at the solution's own time.

| Leg | Setup | What it checks | Measured |
|---|---|---|---|
| P1 | A drop at $T_{sat}$ | $d$ at a quarter, a half and three quarters of $t_{res}$; $T_p = T_{sat}$ | $d$ $2.3\times10^{-9}$, $7.7\times10^{-9}$, $2.4\times10^{-8}$ relative (bound $10^{-6}$); $T_p$ $7\times10^{-12}$ K |
| P2 | P1 past $t_{res}$ | At $t_{res} + \tau_{res}$ and $+ 3\tau_{res}$: $d = d_{res}$, $T_p$ on the residue's relaxation | $d$ 0.023 % below $d_{res}$ (bound 0.1 %); $T_p$ 0.29 and 0.010 K (bounds 1.9 and 0.25 K) |
| P3 | A drop at 2000 K | Ignition at 4.896 ms, then $d$ and $T_p$ at $t_{ign} + t_{res}/2$ | $T_p$ 0.44 K above $T_{sat}$ (bound 1.3 K); $d$ $1.6\times10^{-4}$ relative (bound $4.2\times10^{-4}$) |
| P4 | P1's midpoint with IG and AG | The same numbers | $7.7\times10^{-9}$, as P1 |
| P5 | `gas-phase = gas`, and six inputs that must stop | `INPUT/gas-phase.txt` is read; the refusals name their cause | the drop burns as P1; all six refused |
| P6 | The `[ICE-Physics]` keys on a bare `A 1` line; a `d-res` token over the key; a `combustion=none` token | The keys alone burn as P1, the token wins, `none` leaves the drop inert; `d-res = 0` is refused | $7.7\times10^{-9}$, as P1; the drop inert and heating |

In every leg the eight cells agree to $10^{-10}$, and the velocity and the number density do not move.

Notes on the legs:

- **P2 and P3.** The switches act at the end of a stage, which sets their bounds: the residue can end up to
  $1.2\,K\Delta t/(2 d_{res})$ below $d_{res}$, and an ignited drop holds the temperature it reached in that stage,
  up to $\Delta t\,\mathrm{d}T/\mathrm{d}t$ above $T_{sat}$.
- **P5.** The six refusals: two gas species, no gas phase file, two materials, evaporation on the burning
  material, no `d-res`, and no `boiling-temperature`.

## The cell source

`CombustionUnit` calls the cell source on the same aluminium, with the drop 5 m/s slower than the gas ($Re$ = 60)
and `heat-transfer = Ranz-Marshall`.

| Check | Measured |
|---|---|
| The burn rate against the paper's eq. 4, $2\pi(\mu/Pr)\,n\,d\ln(1 + B)(1 + 0.3\,Re^{1/2}Pr^{1/3})$, typed in the test | $9.5\times10^{-16}$ relative |
| The drop's internal energy does not change while it burns, with emissivity 0 and 0.5 | $6\times10^{-17}$ of $n\dot m(c_pT + L_v)$ |
| The heat to the gas is $n\dot m Q$ | exact |
| Below $T_{sat}$, just below $d_{res}$, and in a gas where $B \le 0$: the source of the material without combustion | bit for bit, no heat to the gas |
