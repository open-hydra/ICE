# Multi-block and MPI

Verification of block connections (ATLAS code `101`) and of the MPI parallelisation. The
[crossing-jets](crossing-jets.md) case with the IG closure, which runs without the shock
detector, is split into four blocks, and the result is compared with the single block and
between different numbers of MPI ranks.

## Setup

The $100 \times 100$ mesh is cut at $x = 0.5$ and $y = 0$ into four $50 \times 50$ blocks.
ATLAS BCB finds the matching faces and writes them as connections. The outer faces keep the
single block's boundary conditions, with the two jet inlets on face 1 of the left-hand
blocks.

| Block | Extent in $x$ | Extent in $y$ |
|---|---|---|
| 1 | $[0, 0.5]$ | $[-0.5, 0]$ |
| 2 | $[0.5, 1]$ | $[-0.5, 0]$ |
| 3 | $[0, 0.5]$ | $[0, 0.5]$ |
| 4 | $[0.5, 1]$ | $[0, 0.5]$ |

This is the regression case `test/Doisneau/IG-split4`. Its `ATLAS/` folder holds the BCB
input for the connections.

## Results

| Check | Result |
|---|---|
| Four blocks vs single block, max difference relative to each variable's maximum | $4 \times 10^{-13}$ (round-off) |
| 2, 3 and 4 MPI ranks vs 1 rank | Bit-identical output |
| [Chimera case](chimera.md) (two blocks) on 2 ranks vs 1 rank | Bit-identical output |
| Single-block cases on 2 ranks (one rank idle) vs serial | Bit-identical output |

A connection passes the two cells behind the interface to the neighbour's ghost layers, so
the four-block run uses the same stencil as the single block and matches it to round-off.
The shock detector would break this: a face on a block boundary takes the detector value of
its own block's cell, while inside a block every face takes the one of the cell on its
low-index side, so wherever a jet edge crosses a block boundary the two runs switch that
face differently; with the detector on, the four-block run differs from the single block by
$2.3 \times 10^{-2}$ of the peak density at 5000 iterations. The MPI runs exchange exactly
the values a serial run reads, so the number of ranks does not change the result at all;
the rank rows above and the timings below were measured with the detector on.

Wall-clock time for 5000 iterations, 2 OpenMP threads per rank:

| MPI ranks | Blocks per rank | Time [s] | Speed-up |
|---|---|---|---|
| 1 | 4 | 41.3 | 1.0 |
| 2 | 2 | 20.3 | 2.0 |
| 3 | 2, 1, 1 | 20.4 | 2.0 |
| 4 | 1 | 13.0 | 3.2 |

On 3 ranks one rank still holds two blocks and sets the pace, so the time is the same as on 2
ranks. The work is divided by whole blocks, so the number of ranks that pays off is one
that divides the blocks evenly.
