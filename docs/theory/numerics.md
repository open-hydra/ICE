# Spatial Discretization

## Finite-volume framework

ICE uses a **cell-centred finite-volume** method on multi-block structured hexahedral grids. The discrete form of the governing equations for cell $i$ is:

$$
\frac{d \mathbf{U}_i}{d t} = -\frac{1}{\Omega_i} \sum_{f \in \partial \Omega_i} \hat{\mathbf{F}}_f \cdot \hat{\mathbf{n}}_f \, A_f + \mathbf{S}_i
$$

where $\Omega_i$ is the cell volume, $A_f$ the face area, $\hat{\mathbf{n}}_f$ the outward unit normal, and $\hat{\mathbf{F}}_f$ the numerical flux at face $f$.

## Convective flux reconstruction

The numerical flux is computed from left and right reconstructed states at each face using an upwind scheme consistent with the monokinetic (pressureless) system. Flux limiters are applied to suppress oscillations near discontinuities.

## Metric tensors

Grid metrics (face areas, cell volumes, unit normals) are computed once at setup from the physical coordinates and stored. The mapping from computational to physical space is evaluated per face.

## Source-term treatment

Stiff drag and heat-transfer source terms are treated **point-implicitly** to relax the stability constraint on the time step without requiring a full implicit solve.

---
