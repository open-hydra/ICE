# Performance Notes

## Where the time goes

A step is dominated by the flux loop, which is swept once per direction per
Runge-Kutta stage per family. Everything else — sources, ghost fill, residual, state
update — is a single pass over the cells.

So the first-order levers are the ones that change how many flux sweeps happen:

| Lever | Effect |
|---|---|
| `time` | Three stages cost about 1.5× two |
| Number of families | Linear: each is a full independent sweep |
| `space-reconstruction` | First order skips the limiter but not the sweep; the saving is modest |
| `dt-max` | Often the real cost driver. If the ceiling binds rather than the CFL condition, every step is smaller than it needs to be — see [Time Integration](../theory/time-integration.md#the-time-step) |

Check `dt-max` first on a slow run: compare the reported `Delta t` against
`cfl * dt-max`. If they are equal, the ceiling is what is setting the pace, and raising
it costs nothing in stability.

## OpenMP

Threading is over the cells of a block, not over blocks, so a single large block
threads as well as many small ones. Two details shape the scaling:

- The interior flux loop is swept in two passes, odd faces then even, because a face
  accumulates into the cells on both sides of it. The barrier between the passes is
  what makes the result independent of the thread count, and it is also a
  synchronisation point per direction per stage.
- The boundary flux loop is `SCHEDULE(DYNAMIC)` over boundary cells, whose cost varies
  a great deal by type — a chimera cell blends several donors, an extrapolation cell
  copies.

Thread counts beyond a few hundred cells per thread stop paying: the loop bodies are
short and the barriers are frequent.

## MPI

Blocks are assigned whole, largest first, to the least-loaded rank. Two consequences
follow directly:

- **Ranks beyond the block count idle.** The startup line `MPI partition: N blocks over
  R ranks, balance X% of ideal` is the thing to read. A balance well under 100% means
  the blocks are unequal, and the fix is in the mesh, not in the launch.
- **Every rank allocates the whole domain**, and only updates its own blocks. Memory
  per rank therefore does not fall as ranks are added, which caps how many ranks fit on
  a node for a large mesh.

The halo exchange runs once per Runge-Kutta stage and carries only the interior cells
that a remote ghost cell reads. Its size is printed at startup as
`MPI halo: N cells exchanged per ghost fill`. Persistent requests are set up once, so
the per-stage cost is the start/wait pair plus the transfer.

Output is gathered to rank 0, which writes alone. For a large mesh written often, this
is a serial section in an otherwise parallel run — lowering `sol-diter` is expensive in
a way the solve itself is not.

## Profiling

```bash
./install.sh build --compilers=gnu --use-openmp
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build build --parallel
```

`gprof`, `perf` and Intel VTune all work on the result. Note that the reported
"Time of operation" in ICE's own output is CPU time divided by the thread count, not
wall-clock, so it is not a substitute for timing the process.
