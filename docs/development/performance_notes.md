# Performance Notes

This page is about what makes one ICE process fast, and how the parallel code is built.
For choosing a parallel configuration and for the measured scaling, see
[Parallel Execution](../user/parallel.md).

## Where the time goes

A step is dominated by the flux loop, swept once per direction, per Runge-Kutta stage,
per family. Everything else — sources, ghost fill, residual, update — is a single pass
over the cells. The settings that matter most are those that change the number of flux
sweeps:

| Setting | Effect |
|---|---|
| `time-scheme` | Three stages cost about 1.5× two |
| Number of families | Linear: each family is a full sweep of its own |
| `space-reconstruction` | First order skips the limiter, not the sweep: a modest saving |
| `dt-max` | Often the real cost: if the ceiling binds instead of the CFL condition, every step is smaller than it needs to be |

On a slow run, check `dt-max` first: if the reported `Delta t` equals `cfl * dt-max`,
the ceiling sets the pace, and raising it costs nothing in stability
([Time Integration](../theory/time-integration.md#the-time-step)).

## OpenMP

The threads share the cells of a block, not the blocks, so one large block threads as
well as many small ones. Which thread computes a cell never changes the arithmetic:
the results are bit-for-bit the same for any thread count.

**The interior flux loop.** Each thread owns a tile of the block — whole k-planes when
the block has at least as many planes as threads, pieces of a plane otherwise. It
computes every face flux of its tile once, into private buffers, and adds the two faces
of each direction to each cell in a fixed order, so the residual does not depend on the
thread count. The shock-detector weights are a pass of their own over the same tiles.

**The boundary flux loop** runs over the boundary cells, each cell on one thread
accumulating its records in face order (a 3-D block corner carries three). The cost per
cell varies a lot: a chimera cell blends several donors, an extrapolation cell copies.

**Scheduling.** The two kinds of loop are scheduled differently:

- the flux tiles are `SCHEDULE(DYNAMIC, 1)`, with `ICE_TILES_PER_THREAD` tiles per thread
  (default 4): the cloud fills only part of the domain, so equal cell counts are not equal
  work, and one static tile per thread left a rank of 80 threads waiting at its barriers
  for a quarter of the sweep. The deep tiling, with private plane buffers, is used only
  while the buffers fit 16 MiB together; above that, the k-faces go to a shared array;
- the per-cell loops (sources, residual, update, time step) are `SCHEDULE(RUNTIME)`, and
  ICE sets that schedule to static at setup, so `OMP_SCHEDULE` has no effect on it (a
  dynamic schedule there cost up to 25× at 80 threads). `ICE_OMP_SCHEDULE=<kind>[,<chunk>]`
  overrides it for measurements.

Beyond a few hundred cells per thread, more threads stop paying: the loop bodies are
short and the barriers frequent.

### First touch

A page of memory lives on the socket of the thread that first writes it.
`Allocate_Block` therefore writes the block fields from a parallel loop over the
k-planes, so each plane starts on the socket of the thread that will sweep it. Written
by one thread, a whole block would start on one socket: a rank of 80 threads over four
sockets would then stream one socket's memory, and its first minute of iterations would
cost several times their steady value while the kernel moves the pages. Ranks of 20
threads pinned to one socket do not need this.

The `ICE Timing` lines show such a start-up transient window by window; a single
end-of-run number hides it.

### Thread groups

A rank whose threads span several sockets and which holds many blocks shares every
block's pages among all the sockets, and the kernel keeps migrating them.
`Mod_ThreadGroups` cuts the threads into groups — by default one per socket, when the
threads are pinned — and gives each block to one group, balancing the cells
(`lpt_assign`, the same rule as blocks over ranks). Every phase then keeps a block on
its group:

- at setup, a group's threads write its blocks' arrays first, so the pages land on its
  socket;
- the per-cell phases give each thread a slice of its group's blocks;
- the ghost fill and the boundary fluxes take each record on a thread of the block's
  group;
- the flux kernel runs in rounds, round $r$ sweeping the $r$-th block of every group at
  once (`compute_flux_grouped`).

`ICE_THREAD_GROUPS=<n>` forces `n` groups, and `1` turns them off. A rank with fewer
blocks than groups uses one group. With the groups, one rank of 80 threads on a
24-block cut runs about as fast as four ranks of 20 threads on the same cut.

Keep large `select case` loop bodies, such as those of the ghost fill, in procedures of
their own: the compiler vectorises more of a body compiled as its own procedure than of
the same body inlined in an OpenMP region.

### The race probe

`Equiv3DGroups` cannot see a missing barrier between the kernel's k-faces and its sweep:
on its small blocks, a sweep tile never reads a seam before its writer has finished.
That barrier has its own test, on a special build outside the `fast` tier and the CI:

```bash
cmake -B build-race -DICE_RACE_PROBE=ON -DICE_BIN_DIR=$PWD/bin-race && cmake --build build-race
ctest --test-dir build-race -L race      # Equiv3DGroupsRace: ~40 s (ifx, RELEASE), ~6 min (gfortran, TESTING)
```

In that build, every round of `compute_flux_grouped` fills the group's shared k-face
planes with NaN and delays one writer, so a missing barrier makes a thread read a NaN
and stops the run. It covers that barrier only. Run it before merging a change to the
grouped kernel; without the option, the normal code is compiled.

## MPI

The blocks go whole to the ranks, largest first, each to the least-loaded rank. So:

- **Ranks beyond the block count are idle.** Read the startup line
  `MPI partition: N blocks over R ranks, balance X% of ideal`: a balance well under 100 %
  means unequal blocks, and the fix is in the mesh.
- **Every rank holds the whole domain** and updates only its own blocks, so the memory
  per rank does not fall as ranks are added.

**The halo exchange** runs once per family per Runge-Kutta stage, and carries the
interior cells a remote ghost cell of that family reads. The families are independent:
a record of family $p$ reads only $p$'s cells. Each family and grid level has one
schedule and one message tag; a cell read by several records is sent once; persistent
requests are set up once, so each exchange costs a start and a wait, the transfer, and
the packing, which the threads share.

Every collective and message of ICE goes through `ice_comm` (`Mod_MPI`): the world,
unless a host program passes its own communicator to `mpi_init_env`.

One consequence for code that reads a ghost cell outside the stage loop (nothing in the
solver does): at the end of a step, family $p$'s ghosts hold the values of $p$'s last
stage.

**Output** is gathered to rank 0, which writes alone: a serial section in an otherwise
parallel run. A small `sol-diter` costs time the solve itself does not.

## Profiling

```bash
cmake -B build -DCMAKE_BUILD_TYPE=TESTING -DCMAKE_Fortran_FLAGS=-g
cmake --build build --parallel
```

The build types are `RELEASE`, `TESTING` (`-O2`) and `DEBUG` (`-O0`, bounds checks);
`RelWithDebInfo` is refused. `gprof`, `perf` and Intel VTune all work on the result, and
`[ICE-IO] timers = true` prints the wall time of each phase and region of the
iteration.

ICE's own "Time of operation" is CPU time divided by the thread count, not wall-clock
time: it does not replace timing the process.
