# Time Integration

## Runge–Kutta scheme

ICE advances the solution in time using an explicit **multi-stage Runge–Kutta** scheme. The number of stages `nrk` is configurable in `input.ini`.

A generic $s$-stage scheme advances from $t^n$ to $t^{n+1} = t^n + \Delta t$ as:

$$
\mathbf{U}^{(0)} = \mathbf{U}^n, \qquad
\mathbf{U}^{(k)} = \sum_{j=0}^{k-1} \alpha_{kj} \mathbf{U}^{(j)} + \beta_{kj} \Delta t\, \mathcal{R}(\mathbf{U}^{(j)}), \qquad
\mathbf{U}^{n+1} = \mathbf{U}^{(s)}
$$

where $\mathcal{R}$ is the residual operator and the coefficients $\alpha_{kj}$, $\beta_{kj}$ define the scheme.

## Time-step control

### Steady-state mode (local time stepping)

Each cell uses its own maximum stable $\Delta t_i$ based on the local CFL condition:

$$
\Delta t_i = \text{CFL} \cdot \frac{\Omega_i}{\sum_f |\mathbf{u}_g \cdot \hat{\mathbf{n}}_f| A_f}
$$

This accelerates convergence to steady state by allowing larger steps in smooth regions.

### Time-accurate mode (global time stepping)

A single global time step is taken:

$$
\Delta t = \text{CFL} \cdot \min_i \Delta t_i
$$

This ensures temporal accuracy and monotone time advancement.

---
