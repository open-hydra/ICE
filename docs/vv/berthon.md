# Berthon Riemann Problems

Three one-dimensional Riemann problems for the anisotropic Gaussian (AG) closure,
compared against their analytical wave patterns. They are the only cases in the suite
in which the AG system carries genuine waves: the [crossing jets](crossing-jets.md) are
free streaming, and the [code-verification](verification.md) cases use the monokinetic
closure, which carries no waves at all. What they exercise is the hyperbolic core — the Rusanov
flux, the MUSCL reconstruction with the MC limiter, and the AG eigenstructure — on
solutions with shocks, contacts and rarefactions.

The exact solutions are tabulated in `test/Berthon/Results/<case>/<case>_exact.tec` as
piecewise-linear profiles of

$$
\rho_p,\quad u_p,\quad v_p,\quad P_{11},\quad P_{22},\quad
\det P = P_{11}P_{22} - P_{12}^2 ,
$$

the last being the quantity the Gaussian closure carries across the contact.

## Setup

Common to the three cases: a $[-0.5, 0.5]$ m tube of 500 cells, one cell across $y$ and
$z$, no gas (0-way coupled), extrapolation at both ends and nothing at the other four
faces. MUSCL with the `mc` limiter, `RK2` in time, Rusanov, $\mathrm{CFL} = 0.5$.

A fourth case, `SCS-hlle`, is the SCS problem through the HLLE flux instead — the only
case that exercises the `riemann-solver` key. HLLE is the less dissipative of the two and its
density error is lower, $3.2\times10^{-3}$ against Rusanov's $4.3\times10^{-3}$.

| | $\rho_p$ | $u_p$ | $v_p$ | $P_{11}$ | $P_{12}$ | $P_{22}$ | $t_{\text{end}}$ |
|---|---|---|---|---|---|---|---|
| **SCS** left  | 1 | 1 | 1 | 1 | 0 | 1 | 0.125 s |
| **SCS** right | 1 | −1 | −1 | 1 | 0 | 1 | |
| **RCS** left  | 1 | 0 | 0 | 2 | 0.05 | 0.6 | 0.125 s |
| **RCS** right | 0.125 | 0 | 0 | 0.2 | 0.1 | 0.2 | |
| **RCR** left  | 2 | −0.5 | −0.5 | 1.5 | 0.5 | 1.5 | 0.150 s |
| **RCR** right | 1 | 1 | 1 | 1 | 0 | 1 | |

$P_{33} = 10^{-6}$, $P_{13} = P_{23} = 0$, $w_p = 0$ and $T_p = 300$ K everywhere; the
number density follows from $\rho_p$ at $\rho_{al} = 2700$ kg/m³ and $R_p = 10^{-6}$ m.
The initial and boundary conditions are written by each case's `IBCB.f90`.

## SCS — shock, contact, shock

Two states colliding head-on along the diagonal. The compression throws a shock each
way, with a contact left between them.

<figure>
  {% include "vv/images/berthon-scs.svg" %}
</figure>

Both shocks sit at the right place and are captured in two to three cells. The visible
defect is the spike at $x = 0$ in $P_{22}$ and $\det P$: a start-up error laid down at
the initial discontinuity, which the contact then carries without spreading. It is the
largest single contribution to this case's error and the reason its tolerance is looser
than the other two.

## RCS — rarefaction, contact, shock

The Gaussian counterpart of a Sod problem: a fourfold pressure ratio expanding to the
right, with unequal shear on the two sides.

<figure>
  {% include "vv/images/berthon-rcs.svg" %}
</figure>

## RCR — two rarefactions

Two states flying apart, which opens a fan on each side of the contact. $P_{22}$ keeps
the plateau structure that the shear $P_{12}$ produces across the contact.

<figure>
  {% include "vv/images/berthon-rcr.svg" %}
</figure>

## What is measured

The solutions are discontinuous, so a pointwise comparison would only report how many
cells the scheme smears a shock over. Each case's `verify.py` therefore takes the $L_1$
error over the domain, normalised by the range of the exact profile:

| $L_1$ error | SCS | RCS | RCR |
|---|---|---|---|
| $\rho_p$ | $4.3\times10^{-3}$ | $2.7\times10^{-3}$ | $1.9\times10^{-3}$ |
| $u_p$ | $1.5\times10^{-3}$ | $1.9\times10^{-3}$ | $7.4\times10^{-4}$ |
| $v_p$ | $4.4\times10^{-3}$ | $7.8\times10^{-3}$ | $3.4\times10^{-3}$ |
| $P_{11}$ | $2.8\times10^{-3}$ | $2.1\times10^{-3}$ | $3.5\times10^{-3}$ |
| $P_{22}$ | $1.5\times10^{-2}$ | $5.3\times10^{-3}$ | $7.9\times10^{-3}$ |
| $\det P$ | $1.6\times10^{-2}$ | $2.3\times10^{-3}$ | $2.8\times10^{-3}$ |
| **tolerance** | $2.5\times10^{-2}$ | $1.5\times10^{-2}$ | $1.5\times10^{-2}$ |

The tolerances are set at roughly twice the measured error: loose enough that a change
in how many cells a shock spreads over will not trip them, tight enough that a wrong
wave speed, a missing wave or a wrong intermediate state — all of which are $O(10^{-1})$
on this scale — will.

!!! note "The reference is a polyline"
    The tabulated exact solution stores a rarefaction as a straight segment between the
    two edges of the fan, while the real fan is curved — visible in $P_{11}$ and
    $\det P$ in the RCR figure, where ICE bows away from the grey line inside the fan.
    Part of the error reported for RCS and RCR is therefore the reference's own
    coarseness, not the solver's, and the numbers above are upper bounds.

!!! warning "The step is set by `dt-max`, not by the CFL number"
    At $\mathrm{CFL} = 0.5$ on this mesh the stability limit is $\Delta t \approx
    1.8\times10^{-4}$ s at the start, but the default `dt-max` holds the step at
    $10^{-4}$ s, so each case takes almost twice as many steps as it needs.
    They still run in about a second each. Raising `dt-max` will change these numbers
    slightly.

## Running them

```bash
ctest --test-dir build -R Berthon --output-on-failure
```

The three Rusanov cases also carry the `fast` label, so the pre-push hook runs them: they are the only
fast-tier cases that touch the AG closure. The figures on this page come from the same
runs:

```bash
cd test/Berthon
python3 plot_vv.py                  # writes docs/vv/images/berthon-*.svg
```
