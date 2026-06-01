# Initial Conditions

ICE reads initial conditions from files in the `ic/` directory of the case. One file is required per block per particle group.

## File format

Initial condition files follow the same binary/ASCII format as ORION output. Each file contains the conservative-variable vector $\mathbf{U}$ at every cell of the block:

$$
\mathbf{U} = (\alpha_p \rho_p,\; \alpha_p \rho_p u,\; \alpha_p \rho_p v,\; \alpha_p \rho_p w,\; \alpha_p \rho_p e_p,\; \alpha_p)
$$

## Naming convention

Files are named following the convention defined in `input.ini`. By default:

```
ic/block_<B>_group_<G>.dat
```

where `<B>` is the 1-based block index and `<G>` is the particle group index.

---
