# Performance Notes

## Profiling

Build with debug symbols and use `gprof` or Intel VTune:

```bash
./install.sh build --compilers=gnu --build-type=RelWithDebInfo
gprof bin/ICE gmon.out > profile.txt
```

## Known bottlenecks

- **Ghost-cell exchange** (MPI runs): the `ICE_Mod_GhostExchange` module performs persistent MPI communications at every Runge–Kutta stage. Reducing the number of RK stages (`nrk`) can significantly cut communication overhead for large rank counts.
- **Source-term evaluation**: drag and heat-transfer loops over all cells and all particle groups. For cases with many groups (`ngroups > 10`) this dominates wall-clock time.

## OpenMP scaling

OpenMP parallelism is applied at the block level. For best scaling, ensure the number of threads does not exceed the number of cells per block divided by the minimum recommended chunk size (~1000 cells/thread).

---
