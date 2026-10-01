# Performance Notes

This page is about what makes a single ICE process fast. For choosing a parallel
configuration and for the measured scaling curves, see
[Parallel Execution](../user/parallel.md).

## Where the time goes

A step is dominated by the flux loop, which is swept once per direction per
Runge-Kutta stage per family. Everything else — sources, ghost fill, residual, state
update — is a single pass over the cells.

So the first-order levers are the ones that change how many flux sweeps happen:

| Lever | Effect |
|---|---|
| `time-scheme` | Three stages cost about 1.5× two |
| Number of families | Linear: each is a full independent sweep |
| `space-reconstruction` | First order skips the limiter but not the sweep; the saving is modest |
| `dt-max` | Often the real cost driver. If the ceiling binds rather than the CFL condition, every step is smaller than it needs to be — see [Time Integration](../theory/time-integration.md#the-time-step) |

Check `dt-max` first on a slow run: compare the reported `Delta t` against
`cfl * dt-max`. If they are equal, the ceiling is what is setting the pace, and raising
it costs nothing in stability.

## OpenMP

Threading is over the cells of a block, not over blocks, so a single large block
threads as well as many small ones. Two details shape the scaling:

- The interior flux loop gathers: each thread owns a tile of the block (a range of
  whole k-planes when the block has at least as many planes as threads, pieces of a
  plane otherwise), computes every face flux of its tile once into private line, row
  and plane buffers, and applies the two faces of each direction to each cell in the
  order the earlier odd-then-even scatter passes gave them — so the residual is bit
  for bit the scatter form's, and independent of the thread count. Per block and stage
  there are three worksharing loops (the shock-detector pass, the faces on the tile
  seams, the sweep) where the scatter form had a pass per parity per direction.
- The boundary flux loop is `SCHEDULE(DYNAMIC, 16)` over the rank's boundary cells of
  the family, each cell accumulating its records in table order on one thread — a
  block corner in 3-D carries three — so 2-D and 3-D runs are bit-identical across
  thread counts (`test/fast/equiv3d`: RED on the record-parallel `ATOMIC` loop, GREEN
  since 2026-09-30). The cost per cell varies a great deal by type: a chimera cell
  blends several donors, an extrapolation cell copies.

Thread counts beyond a few hundred cells per thread stop paying: the loop bodies are
short and the barriers are frequent.

The two loop classes of the step are scheduled differently, and neither choice touches a
result (which thread sweeps a cell is not arithmetic):

- the **flux kernel's tile loops** are `SCHEDULE(DYNAMIC, 1)` over `ICE_TILES_PER_THREAD`
  tiles per thread (default 4), because the particle cloud fills only part of the domain
  and equal cell counts are not equal work: with one static tile per thread a rank of 80
  threads spent 24 % of the sweep's CPU time at its barrier (VTune, monolith, 2026-10-01).
  The kernel takes its deep tiling (private plane buffers) only while the threads' buffers
  together fit `plane_buffer_budget` (16 MiB); above it, the shallow tiling computes every
  k-face once into the shared array.
- the **per-cell loops** (`zero_residual`, `compute_source`, `compute_residual`,
  `state_update`, `state_copy`, `set_dt_global`, `compute_dt`) are `SCHEDULE(RUNTIME)`, and
  ICE sets that schedule itself at set-up: static. `OMP_SCHEDULE` therefore has no effect on
  ICE (a dynamic schedule on these loops dispatches per cell and cost up to 25× at 80
  threads); `ICE_OMP_SCHEDULE=<kind>[,<chunk>]` is the measurement override.

On the 192³ case at one rank of 80 threads (monolith, job 303746, steady state over the last
30 of 60 iterations, two runs each) the two together run 0.303 s per iteration against 0.343
for the kernel before them: −12 %.

**First touch.** `Allocate_Block` writes the block fields (`prim`, `prim_old`, `tau`, `beta`,
`source`, `residual`) from a static parallel loop over the k-planes, so each plane's pages
start on the socket of the thread that sweeps it. Written serially, a whole block starts on
the master's socket: a rank of 80 threads over four sockets then streams one socket's memory
(the whole footprint bound to one node runs 4.15 s per iteration instead of 0.35, job 303745),
and until the kernel's NUMA balancer has moved ~1.7 GB — 45 to 60 s — the iteration costs
6–10× its steady value; on a busy node it may never recover (0.38 s instead of 0.30). With
the parallel first touch every run is flat from its first iteration (jobs 303746). Ranks of
20 threads pinned to one socket do not need it, which is why the multi-rank layouts never
showed the effect. The windows of the `ICE Timing` lines show it; a last-window number does
not.

**Many small blocks on one rank** remain open. On the 24-block cut of the 192³ case (24 × 64 ×
192 blocks) at one rank of 80 threads the current tree runs 0.690 s per iteration against 0.640
for the kernel before it (job 303748): with planes far smaller than a 2 MB huge page, the
plane-by-plane first touch leaves pages shared by threads of two sockets, and the NUMA balancer
keeps migrating them (195 000 hinting faults per run against 5 000–8 000 on one block). With the
balancer off and the placement kept (`numactl --membind=0-3`) the same binary runs 0.584 s; the
best placement the balancer itself reached in one run gave 0.44 s (job 303746). Ranks of 20
threads on their own socket, the recommended layout, are not affected.

## MPI

Blocks are assigned whole, largest first, to the least-loaded rank. Two consequences
follow directly:

- **Ranks beyond the block count idle.** The startup line `MPI partition: N blocks over
  R ranks, balance X% of ideal` is the thing to read. A balance well under 100% means
  the blocks are unequal, and the fix is in the mesh, not in the launch.
- **Every rank allocates the whole domain**, and only updates its own blocks. Memory
  per rank therefore does not fall as ranks are added, which caps how many ranks fit on
  a node for a large mesh.

The halo exchange runs once per family per Runge-Kutta stage and carries the interior
cells that a remote ghost cell of that family reads; the ghost fill that follows it
refills that family's ghosts only (a record of family p reads p's cells alone, so the
families are independent). The size printed at startup, `MPI halo: N cells exchanged
per ghost fill`, is the sum over the families — what the set-up fill exchanges; a stage
of family p moves p's share of it. One schedule and one message tag per family and
grid level; a cell that several records read (a block edge seen from two faces, a
chimera donor of several receivers) is listed once per message — the set-up prints how
many records were merged, when any were; persistent requests are set up once, so the
per-fill cost is the start/wait pair plus the transfer, and the packing and unpacking,
which the threads share (every cell owns a slice of the buffer). Every collective and
message of ICE goes through `ice_comm`
(`Mod_MPI`), the world unless a host program passes its own to `mpi_init_env`. A
consequence for anything that reads a ghost cell outside the stage loop (nothing in the
solver does; a probe placed on a ghost row would): at the end of a step family p's
ghosts hold the values of p's last stage, no longer refreshed by the later families'
stages.

Output is gathered to rank 0, which writes alone. For a large mesh written often, this
is a serial section in an otherwise parallel run — lowering `sol-diter` is expensive in
a way the solve itself is not.

## Profiling

```bash
./install.sh build --compilers=gnu --use-openmp
cmake -B build -DCMAKE_BUILD_TYPE=TESTING -DCMAKE_Fortran_FLAGS=-g && cmake --build build --parallel
```

The build types are `RELEASE`, `TESTING` (`-O2`) and `DEBUG` (`-O0`, bounds checks);
`RelWithDebInfo` is refused by `cmake/SetFortranFlags.cmake`. `gprof`, `perf` and
Intel VTune all work on the result, and `[ICE-IO] timers = true` prints the per-phase
and per-region wall times of the iteration. Note that the reported
"Time of operation" in ICE's own output is CPU time divided by the thread count, not
wall-clock, so it is not a substitute for timing the process.
