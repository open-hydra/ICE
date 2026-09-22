# Chimera Overset

Verification of the chimera (overset) boundary condition, ATLAS code `102`. ATLAS BCB finds
the donor cells that overlap each of the two ghost layers of a chimera face and writes
their volume weights. ICE fills each ghost cell with the weighted blend of its donors,
done in conservative variables, as MOSE does (see
[Boundary Conditions](../user/boundary-conditions.md#chimera)).

All checks use the [crossing-jets](crossing-jets.md) case with the IG closure, and compare
against the same case solved on the original single $100 \times 100$ block.

## Split block with matching cells

The single block is split at $x = 0.5$ into two $50 \times 100$ blocks, and the interface is
declared as chimera with one donor of weight 1 per ghost cell (the cell a block connection
would copy). Each ghost cell then receives exactly the value of one interior cell, so the
chimera path must reproduce the single-block solution.

| Comparison | Max relative difference, all variables |
|---|---|
| Split with chimera vs single block | $1 \times 10^{-11}$ (round-off) |

## Overlapping blocks with non-matching cells

| | Block 1 | Block 2 |
|---|---|---|
| Extent in $x$ | $[0, 0.5]$ | $[0.4, 1.0]$ |
| Cells | $50 \times 100$ | $45 \times 73$ |
| Cell size $\Delta x \times \Delta y$ | $0.0100 \times 0.0100$ | $0.0133 \times 0.0137$ |
| Chimera face | face 2 ($x = 0.5$) | face 1 ($x = 0.4$) |

The blocks overlap by 0.1 m and no cell faces line up, in either direction. Block 1 keeps
the original inlets. ATLAS BCB assigns 1 to 6 donors to each of the 346 chimera ghost
cells, with weights summing to 1, and `chimera.log` reports full coverage for every one of
them.

This is the regression case `test/Doisneau/IG-chimera`. Its `ATLAS/` folder holds the BCB
input used for the donor search.

<figure>
  {% include "vv/images/chimera-overset.svg" %}
</figure>

| Check | Result |
|---|---|
| Block 1 vs single block (mean, relative to peak density) | $6 \times 10^{-5}$ |
| Block 1 vs block 2 inside the overlap (mean / max, relative to peak density) | $6 \times 10^{-4}$ / $3.4 \times 10^{-3}$ |
| Mass flux $\int \rho u\,dy$, block 1 vs single block | $0.00\%$ to $+0.15\%$ |
| Mass flux $\int \rho u\,dy$, block 2 vs single block | $+0.03\%$ to $+0.46\%$ |

The two blocks agree with each other across the overlap, and the jet passes the chimera
interface without a visible step in density or mass flux.

## Matched-resolution control

Block 2 is coarser than the single block, so some of the block 2 deviation above comes from
its grid rather than from the chimera transfer. To separate the two, block 2 was replaced
by a $60 \times 100$ block with the single block's spacing, shifted by half a cell in $x$
($x \in [0.405, 1.005]$). The interface is still non-matching (each ghost cell takes two
donors of weight 0.5), but the resolution is the same on both sides.

| Check | Result |
|---|---|
| Mass flux $\int \rho u\,dy$, both blocks vs single block | within $\pm 0.004\%$ |

So the chimera transfer itself adds a negligible error here, and the 0.03–0.46% seen with
the coarser block 2 comes from its resolution.

!!! note
    The mass flux is computed from cell-centre values of $\rho u$. For the IG closure this
    is not exactly the flux the scheme conserves, so it varies with $x$ even on the single
    block (right panel above). It is used here only to compare runs of the same closure,
    where that effect is common to both.
