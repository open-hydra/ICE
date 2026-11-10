# Parallel Execution

ICE runs in parallel with MPI, OpenMP, or both: MPI shares the mesh blocks among
processes (ranks), and OpenMP shares the cells of a block among the threads of each
rank. The same cores can therefore be used in several ways. On 80 cores:

```text
 1 rank  × 80 threads
 4 ranks × 20 threads
 8 ranks × 10 threads
16 ranks ×  5 threads
```

This page explains how to choose among them, how to launch a run, and how well ICE
scales.

## How ICE divides the work

**MPI shares whole blocks.** The blocks go to the ranks largest first, each to the
rank with the least work so far, and a block is never split. Ranks beyond the number
of blocks have nothing to do, so cut the mesh into at least as many blocks as ranks —
ideally a multiple of them. ICE prints how even the split is:

```text
  MPI partition: 12 blocks over 12 ranks, balance 100.0% of ideal
```

A balance well under 100 % means the blocks are unequal: the fix is in the mesh, not
in the launch command. (ATLAS MDB cuts a mesh into balanced blocks.)

**OpenMP shares the cells of a block.** One large block threads as well as many small
ones, so there is no need to add blocks for the threads. When the threads of one rank
span several sockets and the rank holds several blocks, ICE gives each block to the
threads of one socket (*thread groups*), so a block's memory stays on its socket. This
is automatic when the threads are pinned (`OMP_PROC_BIND`); `ICE_THREAD_GROUPS=<n>`
forces `n` groups, and `1` turns them off.

**Every rank holds the whole domain**, and updates only its own blocks. Before each
Runge-Kutta stage, the ranks exchange the cells their neighbours' ghost cells need;
rank 0 gathers the solution and writes the files. The memory per rank therefore does
not fall as ranks are added: the node's memory, not its cores, limits the number of
ranks ([Memory](#memory)).

**The answer does not depend on the configuration**: a run on any number of ranks and
threads is bit-for-bit identical to a serial one. See
[Multi-block and MPI](../vv/multiblock-mpi.md).

## Choosing a configuration

| Symbol | What it is | How to find it |
| --- | --- | --- |
| `S` | Sockets (NUMA domains) per node | `numactl --hardware` |
| `C` | Cores per socket | `lscpu` |
| `B` | Mesh blocks | ICE's `MPI partition` line, or the ATLAS MDB log |
| `M` | Memory of one rank | the memory of a one-rank run |

For `R` ranks:

1. **Cover every socket**: use at least one rank per socket, or spread the threads of
   one rank over all of them.
2. **Keep each rank's threads on one socket** when the number of ranks allows it.
3. **Check `B ≥ R`**: extra ranks do nothing.
4. **Check `R × M` against the node's memory**: this is usually the real limit.

On a node of 4 sockets × 20 cores:

| Cores | Ranks × threads | Note |
| ---: | ---: | --- |
| 20 | `4 × 5` | one rank per socket |
| 40 | `4 × 10` | one rank per socket |
| 80 | `4 × 20` | one rank per socket: the default choice |
| 80 | `1 × 80` | the least memory; pin the threads (`OMP_PROC_BIND=close`). With several blocks, thread groups make it almost as fast as `4 × 20` |
| 160 | `8 × 20` | two nodes, the same layout on each |
| 240 | `12 × 20` | three nodes, the same layout on each |

Avoid:

- one rank's threads packed on **one socket** while the others sit idle;
- more ranks than blocks;
- more ranks than the node's memory holds;
- interleaving the memory pages over the sockets (`numactl --interleave=all`): it
  makes a run slower, not faster.

## Launching

### One node

With Intel MPI, four ranks of twenty threads, one rank per socket:

```bash
export OMP_NUM_THREADS=20
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n 4 -ppn 4 /path/to/bin/ICE
```

With Open MPI:

```bash
mpirun --map-by socket:PE=$OMP_NUM_THREADS --bind-to core \
       -np 4 /path/to/bin/ICE
```

With OpenMP only, spread the threads over the sockets explicitly — by default they
fill one socket at a time:

```bash
export OMP_NUM_THREADS=80
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
/path/to/bin/ICE
```

### Several nodes

Keep the same layout on every node:

```bash
export OMP_NUM_THREADS=20
mpirun -n 12 -ppn 4 /path/to/bin/ICE     # 3 nodes × 4 ranks × 20 threads
```

## Checking the placement

A run whose threads are not where you meant is still a valid run, only slower, and
nothing warns you. Check before trusting a timing:

```bash
I_MPI_DEBUG=4 mpirun -n 4 -ppn 4 ./bin/ICE 2>&1 | grep -i 'rank.*pin\|domain'
```

Each rank should report its own set of cores, and the sets should cover the sockets
you intended (`numactl --hardware` gives the map of sockets and cores).

## Measuring performance

ICE times its own iterations. Switch the timers on in `[ICE-IO]`; `shell-diter` sets
how often they report:

```ini
[ICE-IO]
timers      = true
shell-diter = 10
```

At 80 cores (`4 × 20`) on the benchmark below, a report looks like this:

```text
 ICE Timing | Iter 40 | 10 iters | wall/iter  3.3361E-01 s | rank min  3.3359E-01 avg  3.3361E-01
            | imbalance  0.0 % | exchange wait    5.9 % | collective wait    2.2 %
 ICE Phases | source  3.7526E-02 s ( 11.2 %) | flux  2.4066E-01 s ( 72.1 %)
            | halo  4.6516E-02 s ( 13.9 %)
 ICE Ranks  | compute/iter max  3.2143E-01 s (rank 0) | min  2.9485E-01 s (rank 2) | mean  3.0659E-01 s | spread   4.8 %
 ICE Cycles | core-cycles/iter  6.9599E+10 | aggregate    208.62 GHz
```

For numbers you can compare:

- run enough iterations to have several reports, and read the **last** one: the first
  includes the start-up;
- switch the solution output off (a large `sol-diter`) and keep `res-diter` large;
- repeat each configuration and keep the fastest run.

### Core-cycles, not only seconds

`core-cycles/iter` is the work of one iteration summed over every rank and thread.
With perfect parallelism it does not change as cores are added, so its ratio,
$C(1)/C(N)$, is the parallel efficiency without the effect of the processor's clock.
That effect is large: on the benchmark node, one busy core runs at 3.87 GHz and
eighty at 2.61 GHz, so a full node loses a third of its speed per core to the chip's
power management, not to ICE.

The counter needs a C compiler at configure time (`ICE_CYCLES`, on by default).
Without it the `ICE Cycles` line is missing and the timers report seconds only.

## How ICE scales

The benchmark is a cube of 192³ = 7,077,888 cells: one family with the IG closure,
one-way coupled to a frozen Taylor–Green vortex, RK2 with MUSCL and the van Leer
limiter, CFL 0.5. The machine has 4 sockets × 20 cores per node, with no
hyperthreading. Each rank count runs its own cut of the mesh.

| Cores | Nodes | Ranks × threads | Time per iteration | Speed-up | Core-cycle efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | `1 × 1` | 13.13 s | 1.0 | 100.0 % |
| 4 | 1 | `4 × 1` | 3.335 s | 3.9 | 98.4 % |
| 8 | 1 | `4 × 2` | 1.726 s | 7.6 | 95.9 % |
| 16 | 1 | `4 × 4` | 0.905 s | 14.5 | 96.2 % |
| 20 | 1 | `4 × 5` | 0.742 s | 17.7 | 95.9 % |
| 40 | 1 | `4 × 10` | 0.428 s | 30.7 | 89.4 % |
| 80 | 1 | `4 × 20` | 0.295 s | 44.5 | 82.1 % |
| 160 | 2 | `8 × 20` | 0.148 s | 89.0 | 82.4 % |
| 240 | 3 | `12 × 20` | 0.101 s | 130.5 | 80.6 % |

The core-cycle efficiency stays near 96 % up to a quarter of a node, is 82 % on a full
node, and holds at 81 % on three nodes. The speed-up is lower than the core count
mostly because of the clock: in seconds, a full node reaches 56 % of linear.

Part of the remaining loss is not parallel overhead: cutting the mesh adds ghost cells,
4.2 % more cells at 4 ranks and 10.4 % at 12, which the one-rank run does not compute.

<figure>
  --8<-- "docs/user/images/scaling-envelope.svg"
<figcaption>The same curve, as measured in September 2026, before the optimisations of October: they made every multi-rank point 10 to 14 % faster.</figcaption>
</figure>

### MPI, OpenMP or both

All three work, and their differences are small compared with the effect of thread
placement. Core-cycle efficiency, measured in September 2026:

| Cores | OpenMP `1 × N` | MPI `N × 1` | Hybrid `4 × N/4` |
| ---: | ---: | ---: | ---: |
| 4 | 99.0 % | 97.5 % | 97.5 % |
| 8 | 96.8 % | 94.7 % | 95.7 % |
| 16 | 95.2 % | 91.3 % | 95.1 % |
| 20 | 94.6 % | — | 95.5 % |
| 24 | — | 86.8 % | — |
| 40 | 90.6 % | — | 89.2 % |
| 80 | 83.4 % | — | 80.8 % |

- **OpenMP alone reads highest because it computes no ghost cells**: it runs the uncut
  mesh. Credited for its ghost cells, the hybrid `4 × 20` is the more efficient.
- **Pure MPI stops early because of memory, not efficiency**: 24 ranks already hold
  143 GB on one node.

<figure>
  --8<-- "docs/user/images/scaling-modes.svg"
<figcaption>The three modes, September 2026.</figcaption>
</figure>

### The rank/thread split

At 80 cores on one node, every split is within 3 % of the fastest, while its memory
varies fifteenfold (September 2026):

| Layout | Ranks per socket | Time per iteration | Memory of the node |
| --- | ---: | ---: | ---: |
| `1 × 80` | — | 0.324 s | 6.2 GB |
| `8 × 10` | 2 | 0.327 s | 45.2 GB |
| `16 × 5` | 4 | 0.333 s | 92.8 GB |
| `4 × 20` | 1 | 0.334 s | 24.5 GB |

**Do not over-tune this.** Choose the number of ranks from the memory and the blocks,
and move on. After the October optimisations `4 × 20` takes 0.293 s and `1 × 80` (one
block) 0.303 s: one rank per socket stays the default.

<figure>
  --8<-- "docs/user/images/scaling-layout.svg"
<figcaption>The rank/thread split at 80 cores, September 2026.</figcaption>
</figure>

### Across nodes

Adding nodes is nearly free. Compared with one full node:

| | Speed-up | Of linear |
| --- | ---: | ---: |
| 160 cores, 2 nodes | 2.00 × | 100 % |
| 240 cores, 3 nodes | 2.93 × | 98 % |

The computation itself scales almost perfectly across nodes; what grows is the
exchange of ghost cells. If wall-clock time matters most, spreading a case over more
nodes costs little efficiency — only node-hours.

<figure>
  --8<-- "docs/user/images/scaling-nodes.svg"
<figcaption>Scaling across nodes, September 2026.</figcaption>
</figure>

## Memory

Every rank holds the whole domain, so the memory of a node grows with the number of
ranks instead of shrinking. On the benchmark (7 M cells):

| Ranks on the node | 1 | 2 | 4 | 8 | 16 | 24 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Memory (GB) | 6.2 | 12.4 | 24.3 | 45.1 | 93.0 | 143.0 |

That is about 6 GB per rank at this size: a node of 187 GB holds about 30 ranks, so pure
MPI on 80 cores (some 480 GB) is out of reach. When memory is tight, use fewer ranks
and more threads per rank.

Each grid level of `[ICE-Multigrid]` adds a full set of arrays on every rank, about an
eighth of the fine grid's in 3-D and a quarter in 2-D, for the whole run.

<figure>
  --8<-- "docs/user/images/scaling-memory.svg"
<figcaption>Memory against the number of ranks, September 2026.</figcaption>
</figure>

## Known limits

- **Memory per rank**: the whole domain on every rank caps the number of ranks per
  node.
- **Blocks**: ranks beyond the number of blocks are idle, and unequal blocks show in the
  `balance` figure. Cut the mesh for the number of ranks you intend to use.
- **Equal cells are not equal work**: even at 100 % cell balance, the ranks' compute
  times differ by up to 6 %. The cost per cell depends on the solution, and a symmetry
  face costs more than a connection.
- **Many small blocks on one rank over several sockets** are slow without thread
  groups, which need pinned threads and at least as many blocks as sockets.
- **Output is serial**: rank 0 writes alone. Keep `sol-diter` large in big parallel
  runs.
- **One problem size**: these measurements are one 7 M-cell case. Larger meshes, more
  families or more physics will behave differently.

## Quick reference

| Goal | How |
| --- | --- |
| Use every socket | at least one rank per socket, or `OMP_PROC_BIND=spread` |
| Keep memory local | each rank's threads on one socket; no page interleaving |
| Save memory | fewer ranks, more threads |
| Avoid idle ranks | at least as many blocks as ranks |
| Scale to more nodes | the same layout on every node |
| Measure | `timers = true`, the last report, the fastest of several runs |
| Measure efficiency | `core-cycles/iter`, not seconds |
| Check the placement | `I_MPI_DEBUG=4` and the core of each rank |
