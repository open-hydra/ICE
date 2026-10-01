# Parallel Execution and Scalability

ICE supports **hybrid MPI + OpenMP parallel execution**. MPI distributes the
mesh blocks across processes; OpenMP parallelises the work inside each process,
over the cells of a block.

The same core budget can therefore be spent in several ways. On 80 cores:

```text
 1 MPI rank  × 80 OpenMP threads
 4 MPI ranks × 20 OpenMP threads
 8 MPI ranks × 10 OpenMP threads
16 MPI ranks ×  5 OpenMP threads
```

---

## 1. Work division

Three properties of ICE decide what a sensible configuration looks like.

**MPI distributes whole blocks.** Blocks are assigned largest first to the
least-loaded rank, and a block is never split. Ranks beyond the block count sit
idle, so the mesh must be cut into at least as many blocks as there are ranks —
and ideally into a multiple of the rank count, or the balance is poor. ICE
prints what it decided at startup:

```text
  MPI partition: 12 blocks over 12 ranks, balance 100.0% of ideal
  MPI halo: <n> cells exchanged per ghost fill
```

A balance well under 100 % means the blocks are unequal, and the fix is in the
mesh rather than in the launch command.

**OpenMP threads over the cells of a block**, not over blocks. A single large
block therefore threads as well as many small ones, and adding blocks purely to
feed threads is never necessary — but see [§11](#11-known-scalability-limits)
for what happens when a rank holds many small blocks at a high thread count.

**Every rank allocates the whole domain** and updates only its own blocks. Per
rank memory does not fall as ranks are added, so the node's memory — not its
cores — is what caps the rank count. See [§10](#10-memory).

The result does not depend on the parallel configuration: output is bit-for-bit
identical across rank counts. See
[Multi-block and MPI](../vv/multiblock-mpi.md).

---

## 2. Selecting a parallel configuration

| Symbol | Description | How to obtain it |
| --- | --- | --- |
| `S` | Sockets (NUMA domains) per node | `numactl --hardware` |
| `C` | Physical cores per socket | `lscpu` |
| `B` | Number of mesh blocks | ICE's `MPI partition` line, or the ATLAS MDB log |
| `M` | Memory per rank for the case | ICE's resident set at one rank |

For a core budget `N`:

1. **Cover all sockets.** Use at least one MPI rank per socket, or place the
   threads of a single rank across the sockets explicitly.
2. **Keep each OpenMP team within one socket** where the rank count allows it.
3. **Check `B ≥ R`.** Ranks beyond the block count do nothing.
4. **Check `R × M` against the node's memory.** This is usually the binding
   constraint, not the core count.

On a node with **4 sockets × 20 cores**:

| CPU cores | Ranks × threads | Note |
| ---: | ---: | --- |
| 20 | `4 × 5` | one rank per socket |
| 40 | `4 × 10` | one rank per socket |
| 80 | `4 × 20` | one rank per socket, one node |
| 80 | `1 × 80` | lowest memory, needs explicit placement; fastest only on the pre-wave baseline — see [§6](#6-strong-scaling) |
| 160 | `8 × 20` | two nodes, same per-node layout |
| 240 | `12 × 20` | three nodes, same per-node layout |

### Configurations to avoid

* An OpenMP team packed onto **one socket** while other sockets sit idle.
* More ranks than blocks.
* A rank count whose memory exceeds the node.
* Interleaving memory pages across sockets (e.g. `numactl --interleave=all`)
  for a single large-thread-count rank — it makes a thread's private buffers
  slower to reach, not faster. See [§6](#6-strong-scaling).

---

## 3. Running ICE in parallel

### Single node

```bash
export OMP_NUM_THREADS=20
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n 4 -ppn 4 /path/to/bin/ICE
```

This launches 4 ranks × 20 threads = 80 cores, one rank per socket, each rank's
threads contained in its own socket.

With OpenMPI the equivalent is:

```bash
mpirun --map-by socket:PE=$OMP_NUM_THREADS --bind-to core \
       -np 4 /path/to/bin/ICE
```

For a pure-OpenMP run spread over all four sockets, pin explicitly — the default
placement fills one socket at a time:

```bash
export OMP_NUM_THREADS=80
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
/path/to/bin/ICE
```

### Multiple nodes

Keep the per-node layout identical on every node:

```text
 80 cores:  1 node  × 4 ranks × 20 threads
160 cores:  2 nodes × 4 ranks × 20 threads/node
240 cores:  3 nodes × 4 ranks × 20 threads/node
```

```bash
export OMP_NUM_THREADS=20
mpirun -n 12 -ppn 4 /path/to/bin/ICE
```

---

## 4. Verifying CPU placement

A run that is silently unpinned looks like a scaling result and is not one.
Check before trusting any timing:

```bash
I_MPI_DEBUG=4 mpirun -n 4 -ppn 4 ./bin/ICE 2>&1 | grep -i 'rank.*pin\|domain'
```

Each rank should report a distinct set of cores, and those sets should cover the
sockets you intended. `numactl --hardware` gives the socket-to-core map to check
against.

---

## 5. Measuring parallel scaling

ICE has its own per-iteration instrumentation. It is **off by default**; turn
it on in the `[ICE-IO]` section of `input.ini`, where `shell-diter` sets how
often it reports:

```ini
[ICE-IO]
timers      = true
shell-diter = 10
```

Read the **last** window, not the wall clock of the batch job. At 80 cores on
the reference benchmark (job 303350, tag `env-C80-4x20`, third of three
repetitions, iteration 40's window — `runs/303350/logs/env-C80-4x20-r3.log`)
it prints:

```text
 ICE Timing | Iter 40 | 10 iters | wall/iter  3.3361E-01 s | rank min  3.3359E-01 avg  3.3361E-01
            | imbalance  0.0 % | exchange wait    5.9 % | collective wait    2.2 %
 ICE Phases | source  3.7526E-02 s ( 11.2 %) | flux  2.4066E-01 s ( 72.1 %)
            | halo  4.6516E-02 s ( 13.9 %)
 ICE Ranks  | compute/iter max  3.2143E-01 s (rank 0) | min  2.9485E-01 s (rank 2) | mean  3.0659E-01 s | spread   4.8 %
 ICE Cycles | core-cycles/iter  6.9599E+10 | aggregate    208.62 GHz
```

For measurements that mean something:

* run enough iterations that the first window is not scored;
* score the **last** window;
* disable solution output (`sol-diter`) and keep `res-diter` high;
* repeat each configuration and take the minimum.

### Use core-cycles, not seconds

`core-cycles/iter` is the total over every rank and thread. Under perfect
parallelism it stays **constant** as cores are added, so it is an efficiency
measure with the processor's frequency scaling already divided out.

This matters more than it sounds. On the reference node one active core turbos
to 3.87 GHz while eighty sit at 2.61 (job 303350: `env-C1-1x1` vs.
`env-C80-4x20`). That factor of 0.67 is the chip's power management, not ICE,
and it accounts for essentially the whole difference between ICE's wall-clock
efficiency at a full node (54.4 %) and its core-cycle efficiency (80.8 %).

The counter needs a C compiler at configure time (`ICE_CYCLES`, on by default).
Without it the `ICE Cycles` line is absent and the timers report wall-clock
only.

---

## 6. Strong scaling

Strong scaling measures how execution time changes when the **same problem** is
solved on more cores.

The reference benchmark is **7,077,888 cells** (192³), one dispersed-phase
family with the IG closure, one-way coupled to a frozen Taylor–Green carrier,
RK2 with MUSCL/vanleer at CFL 0.5. The reference machine has 4 sockets × 20
physical cores per node and no hyperthreading. The campaign and its harness are
Marco Grossi's (job 302995, first full matrix); the table below is job 303350,
the same harness re-run on `origin/merging` `aeb34a5` — the tree as it stood
before this session's fixes — on nodes wn[05-07]. It is the only run in this
campaign with all nine core counts, so it is the source for every table in this
page unless stated otherwise; what the fixes changed is not yet measured at
that resolution — see "Since the baseline" below.

<figure>
  {% include "user/images/scaling-envelope.svg" %}
</figure>

*(this figure predates job 303350 and the fixes below; it has not been
regenerated for the numbers on this page)*

| Cores | Nodes | MPI × OpenMP | Wall/iter | Core-cycle efficiency | Ghost-credited |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | `1 × 1` | 14.516 s | 100.0 % | 100.0 % |
| 4 | 1 | `4 × 1` | 3.721 s | 97.5 % | 101.6 % |
| 8 | 1 | `4 × 2` | 1.913 s | 95.7 % | 99.7 % |
| 16 | 1 | `4 × 4` | 1.011 s | 95.1 % | 99.1 % |
| 20 | 1 | `4 × 5` | 0.823 s | 95.5 % | 99.5 % |
| 40 | 1 | `4 × 10` | 0.475 s | 89.2 % | 92.9 % |
| 80 | 1 | `4 × 20` | 0.334 s | 80.8 % | 84.2 % |
| 160 | 2 | `8 × 20` | 0.167 s | 80.5 % | 85.4 % |
| 240 | 3 | `12 × 20` | 0.117 s | 76.5 % | 84.5 % |

Both efficiency columns are core-cycle efficiencies, `C(1)/C(N)` ([§5](#5-measuring-parallel-scaling)), not
wall-clock ones: in seconds the full node reaches 54.4 %, and the difference is the
chip's clock, 3.87 GHz on one active core against 2.61 GHz on eighty. Core-cycle
efficiency stays above **95 % through a quarter node** (4–20 cores), reaches
**80.8 % at a full node**, and eases only slightly further to **76.5 % at three
nodes**.

**The two efficiency columns differ because the decomposition is real work.**
Each rank count runs its own cut of the mesh, and cutting adds ghost cells:
4.2 % more cells at 4 ranks, 10.4 % at 12 (a fixed property of this mesh's cut
at each rank count — `analyze.py`'s ghost-cell table, shared with the ATLAS MDB
decomposition analysis — not something that varies between jobs). The raw
column charges ICE for them; the ghost-credited column does not. Quote the raw
column when asking "what does this many ranks cost me", and the ghost-credited
column when asking "how well does the parallelisation itself hold up".

### Since the baseline: the 2026-09-30 optimization wave

The full nine-point matrix on the final tree is job 303485 (tree `461eb4f`,
the final solver code plus one change measured neutral and reverted
afterwards), run on the same three nodes (wn[05-07]) as the baseline above,
so the comparison is clean:

| Cores | Nodes | MPI × OpenMP | 303350 wall/iter | 303485 wall/iter | Change | 303350 core-cycle eff. | 303485 core-cycle eff. |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | `1 × 1` | 14.52 s | 13.13 s | -9.5 % | 100.0 % | 100.0 % |
| 4 | 1 | `4 × 1` | 3.721 s | 3.335 s | -10.4 % | 97.5 % | 98.4 % |
| 8 | 1 | `4 × 2` | 1.913 s | 1.726 s | -9.8 % | 95.7 % | 95.9 % |
| 16 | 1 | `4 × 4` | 1.011 s | 0.9046 s | -10.6 % | 95.1 % | 96.2 % |
| 20 | 1 | `4 × 5` | 0.8233 s | 0.7422 s | -9.9 % | 95.5 % | 95.9 % |
| 40 | 1 | `4 × 10` | 0.4753 s | 0.4278 s | -10.0 % | 89.2 % | 89.4 % |
| 80 | 1 | `4 × 20` | 0.3336 s | 0.2952 s | -11.5 % | 80.8 % | 82.1 % |
| 160 | 2 | `8 × 20` | 0.1671 s | 0.1475 s | -11.7 % | 80.5 % | 82.4 % |
| 240 | 3 | `12 × 20` | 0.1174 s | 0.1006 s | -14.3 % | 76.5 % | 80.6 % |

Core-cycle efficiency is anchor cycles over cycles per iteration, each job
against its own anchor. The multi-rank points gained 9.5–14.5 % in time and
1.3–4.1 points of core-cycle efficiency (the halo work of
[§9](#9-scaling-across-nodes) is where the three-node point gained most);
pure MPI at 24 ranks went from 86.8 % to 90.0 %, the 16 × 5 layout from
81.0 % to 83.3 %, one socket of 20 threads from 94.1 % to 92.9 % (0.837 →
0.767 s). The same matrix, at the four points first measured with the reduced
job 303417 (tree `5575f65`, the same nodes):

| Cores | Layout | 303350 (baseline) | 303417 (wave) | Change |
| ---: | --- | ---: | ---: | ---: |
| 1 | `1 × 1` (anchor) | 14.516 s | 13.129 s | −9.5 % |
| 80 | `4 × 20` | 0.3336 s | 0.2957 s | −11.4 % |
| 240 | `12 × 20` | 0.1174 s | 0.1004 s | −14.5 % |
| 80 | `1 × 80` | 0.3241 s | 0.3368 s | +3.9 % |

Every multi-rank point got faster. The single-rank, 80-thread point got
**slower**, and its cause is **open**: two candidate fixes were built, proven
byte-identical and measured on the same nodes, and neither moved it. Slab
tiles that keep the flux kernel's per-thread plane buffers inside the cache
(job 303439) took 6 % off the sweep but left the leg at 0.338 s, cost 31–40 %
on the 20-thread socket legs through an unbalanced tile count and doubled the
many-blocks leg, so they were reverted; a parallel first touch of the block
fields (job 303482) changed nothing at all (0.3373 s) and was reverted too —
the kernel's automatic NUMA balancing very likely places the pages within the
first iterations already. On the current tree the leg pays about 25 % more in
the per-cell streaming phases and 12 % more in the flux sweep than `4 × 20`
does for the same cells; a memory-access profile of `omp-C80-1x80` is the next
step.

Interleaving memory pages does **not** fix this, and makes it markedly worse on
the current kernel. Job 303418 measured `1 × 80` and `1 × 40` on the wave's
tree (`5575f65`) with the pages left where the serial first touch put them
against the pages explicitly interleaved across sockets: interleaving cost
**3.7×** the wall time at `1 × 80` and **1.28×** at `1 × 40`. On the
pre-wave binary (job 303351, same node), the same comparison cost only
8.0–9.5 % — the interleaving penalty grew with the kernel change, likely
because the hot data are the per-thread plane buffers, which a thread allocates
inside the region and which interleaving spreads across sockets instead of
leaving local (a hypothesis: shrinking those buffers, job 303439, did not move
the plain `1 × 80` time). "One rank per socket" ([§2](#2-selecting-a-parallel-configuration))
remains the right default for this reason.

Until the single-rank leg is understood, prefer a multi-rank layout
(`4 × 20` or more) on one node: it is 11 % faster than the baseline where
`1 × 80` is 4 % slower.

**Decomposition quality** (job 303484, greedy cuts of the campaign against
the halo-objective cuts of ATLAS MDB's search, same binary): at 12 ranks the
two are equal (0.1004 vs. 0.1006 s at `12 × 20`; the halo is a few per cent of
the iteration there), at 16 ranks equal within noise, at 24 ranks the search's
cut is 3.4 % faster in pure MPI (`24 × 1`, 0.656 → 0.633 s) and 4.3 % faster at
`4 × 20` with six blocks per rank (0.320 → 0.306 s, exchange wait 6.3 → 3.3 %),
and 20 ranks — which the greedy cutter could not balance — run at 87 % wall
efficiency on the search's 5×2×2 cut.

---

## 7. MPI, OpenMP and hybrid execution

<figure>
  {% include "user/images/scaling-modes.svg" %}
</figure>

*(this figure predates job 303350; it has not been regenerated for the numbers
on this page)*

Job 303350 (tree `aeb34a5`, nodes wn[05-07]), `omp-*` / `mpi-*` / `env-*` rows:

| Cores | Pure OpenMP `1 × N` | Pure MPI `N × 1` | Hybrid `4 × N/4` |
| ---: | ---: | ---: | ---: |
| 4 | **99.0 %** | 97.5 % | 97.5 % |
| 8 | **96.8 %** | 94.7 % | 95.7 % |
| 16 | **95.2 %** | 91.3 % | 95.1 % |
| 20 | **94.6 %** | — | 95.5 % |
| 24 | — | 86.8 % | — |
| 40 | **90.6 %** | — | 89.2 % |
| 80 | **83.4 %** | — | 80.8 % |

All three models work, and the spread between them is small compared with what
placement does. Pure OpenMP reads highest because it runs an undecomposed mesh
and therefore computes no ghost cells at all; per cell of real physics the
hybrid and MPI runs are the more efficient.

Three things to read carefully:

1. **The pure modes are pinned across all four sockets at every core count**, so
   they are compared on equal hardware rather than on whatever the default
   placement gives. The hybrid points are launched the way a user would, with
   `I_MPI_PIN_DOMAIN=omp` and no mask, and below a full node those domains are
   assigned in core order and need not reach every socket ([§4](#4-verifying-cpu-placement)).
   On this campaign's own measurement, though, an explicit spread mask cost no
   more than the default at 20 cores — `sock-C20-4x5-spread` reads 95.6 % against
   the table's 95.5 %, a 0.1-point difference, inside the campaign's own
   ~1.3 % core-cycle noise floor. That does not mean the effect is zero at
   every core count, only that it was not visible here; only the 80-core
   points fill the node and are free of the question by construction.
2. **Pure OpenMP's advantage is the ghost cells, not the threading.** At 80
   cores it computes 4.2–12.5 % fewer cells than any of the hybrid layouts.
   Credit those cells back and the ranking inverts (84.2 % for `4 × 20` against
   83.4 % for `1 × 80`) — see
   [§8](#8-choosing-the-rankthread-split).
3. **Pure MPI's 24-rank point carries 18.8 % more cells** than the reference —
   the greedy cut's ghost overhead at 24 ranks. Credited for it, its
   efficiency reads **103.2 %**: not "in line with the rest", since every
   hybrid/MPI point up to 16 ranks already credits to 100 % or a little above
   (101.6 %, 100.6 %, 102.8 %), so this is not a 24-rank-specific effect —
   plausibly the smaller per-rank working set fitting cache better than the
   single-block anchor does. What ends the pure-MPI curve well short of 80
   cores is memory ([§10](#10-memory)), not efficiency: 24 ranks already hold
   143 GB on one node.

---

## 8. Choosing the rank/thread split

<figure>
  {% include "user/images/scaling-layout.svg" %}
</figure>

*(this figure predates job 303350; it has not been regenerated for the numbers
on this page)*

At 80 cores on one node, job 303350, `lay-*` rows plus `omp-C80-1x80`:

| Layout | Ranks per socket | Wall/iter | vs best | Node memory |
| --- | ---: | ---: | ---: | ---: |
| `1 × 80` | — | 0.3241 s | best | 6.21 GB |
| `8 × 10` | 2 | 0.3271 s | +0.9 % | 45.16 GB |
| `4 × 20` | 1 | 0.3339 s | +3.0 % | 24.46 GB |
| `16 × 5` | 4 | 0.3330 s | +2.7 % | 92.82 GB |

**Do not over-tune this.** The whole family spans 3.0 %, which is close to the
measurement floor, while the memory it costs spans a factor of fifteen. Choose
the rank count on memory and on the block count, and move on.

This ranking is the pre-wave baseline's. The wave that follows it in
[§6](#6-strong-scaling) changes `1 × 80` specifically — by 3.9 % on its own
tree, and by far more if pages get interleaved — so re-check §6 before reading
`1 × 80` as current advice.

---

## 9. Scaling across nodes

<figure>
  {% include "user/images/scaling-nodes.svg" %}
</figure>

*(this figure predates job 303350; it has not been regenerated for the numbers
on this page)*

Adding nodes is nearly free on the baseline (job 303350). Referred to the
80-core point on one node:

| | Speed-up | Of linear |
| --- | ---: | ---: |
| 160 cores, 2 nodes | 1.996 × | **99.8 %** |
| 240 cores, 3 nodes | 2.843 × | 94.8 % |

The first node boundary costs nothing measurable. Splitting the iteration by
phase locates the remaining ~5 % — using the `ICE Phases` line's `source +
flux` as "compute" and its `halo` term as "halo", best rep of three at each
point: from 80 to 240 cores compute goes from 0.2782 s to 0.0942 s (2.96× of a
3× ideal, **98.5 % of linear**) while the halo phase goes from 0.0465 s to
0.0236 s (1.97× of 3×, **65.6 % of linear**). The halo phase is the limit, not
the interconnect — the Timing line's `exchange wait` is the same 5.9–6.1 % at
80 and 160 cores and only grows at 240 (10.0 %). See
[§11](#11-known-scalability-limits) for how much the wave has since narrowed
the halo phase itself.

The practical consequence: **if wall-clock time is what matters, spreading a
fixed case over more nodes is close to free in efficiency.** It costs node-hours
in proportion, so it is a trade against allocation budget, not against
scalability.

---

## 10. Memory

<figure>
  {% include "user/images/scaling-memory.svg" %}
</figure>

*(this figure predates job 303350; it has not been regenerated for the numbers
on this page)*

Every MPI rank allocates the **whole domain** and updates only its own blocks,
so the node footprint grows linearly with rank count instead of shrinking with
the decomposition. Job 303350, `env-C1-1x1` and the pure-MPI `mpi-*` rows:

| Ranks | 1 | 2 | 4 | 8 | 16 | 24 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GB/node | 6.21 | 12.36 | 24.32 | 45.14 | 92.95 | 143.01 |

That is **6.2 GB for a single rank and ~6.0 GB per rank** at 24 ranks (the
largest single rank measured 6.8 GB), at 7 M cells, whatever the rank count
otherwise. A 187 GB node therefore holds about **30 ranks** at this problem
size, and pure MPI at a full node is impossible: `80 × 1` would need some
480 GB.

This is why the rank count should be chosen on memory. Fewer ranks, more
OpenMP threads per rank is the right lever when memory is tight — and on the
pre-wave baseline it cost nothing in speed to do so ([§8](#8-choosing-the-rankthread-split)).
The optimization wave changes that at the single-rank end: see
[§6](#6-strong-scaling) before assuming a `1 × N` layout is still free.

---

## 11. Known scalability limits

**Memory per rank.** Domain-wide structures live on every rank. Very large rank
counts exhaust memory before they exhaust cores.

**Block count.** Ranks beyond the number of blocks idle, and an uneven split
shows up directly in the `balance` figure at startup. Cut the mesh for the rank
count you intend to use.

**Many small blocks on one rank.** Threading is over cells, not blocks
([§1](#1-work-division)), but the two are not independent at a high thread
count. Cutting the 192³ case into 24 blocks of 24×64×192 instead of the
default 4 blocks of 96×96×192, then running that cut on one rank of 80 threads
(job 303417, tag `blk-C80-1x80-R24`), makes every per-cell phase 2–8× slower
than the 4-block cut on the same rank — `zero` 70.6 ms vs. 8.7 ms, `source` 87
vs. 37, `update` 67 vs. 24, `bound` 21 vs. 3.8, `ghost` 23 vs. 4 — while the
flux sweep itself is essentially unchanged (0.228 s vs. 0.220 s). The `C4b`
NOWAIT change only took `zero` from 107 ms to 70.6 ms, so it is not a barrier
cost, and the per-thread chunk is the same ~3.7k cells (236 kB) in both
layouts. Something in how 80-thread worksharing splits many small blocks costs
roughly 2–3 ms per block per phase; this is **unresolved** — a threading
profile of that leg is the next step. At `4 × 20` the same 24-block cut costs a
much smaller +16 % (job 303350, `blk-C80-4x20-R24` vs. `env-C80-4x20`); at
`1 × 80` it costs +54 %. Practical advice: as few blocks per rank as the
balance allows, and it matters far more at `1 × N` than at a multi-rank split.

**Cell balance is not work balance.** At reported 100 % cell-count balance
at every rank count in this benchmark, ICE still measures up to 6.3 % spread
in per-rank compute time (job 303350, `env-C240-12x20`). Two reasons: the
reconstruction retries on unphysical states, so per-cell cost is data-dependent;
and a rank's mix of physical boundaries and cut connections is not balanced by
balancing cells — a symmetry condition costs far more than a connection's
straight copy.

**The halo phase.** On the pre-wave baseline (job 303350) it is 20.1 % of the
iteration at 240 cores and scales at 65.6 % of linear from 80 to 240 cores
([§9](#9-scaling-across-nodes)) — the limit on multi-node runs, not the
interconnect. The wave's halo work (`C9a`: unique halo cells and MPI statuses;
`C9b`: threaded pack/unpack) cuts this sharply: at `12 × 20`, comparing job
303387 (before, tree `5265964`) to job 303402 (after, tree `e7cdef7` =
C1+C9a+C9b, best rep of two at each), `pack` goes from 5.9 % of the iteration
to 0.7 %, `wait` from 8.1 % to 4.2 %, `unpack` from 4.3 % to 0.7 %, and the
Timing line's `exchange wait` from 11.6 % to 3.8 %. On the C4a+C4b tree (job
303417, same nodes as the baseline) the halo phase itself is down to 10.8 % of
the iteration at 240 cores — about half the baseline's 20.1 %.

**Serial output.** Solution output is gathered through rank 0, which writes
alone. Keep `sol-diter` high for large parallel runs.

**Grid levels cost memory.** Each level in `[ICE-Multigrid]` adds a complete set
of block arrays and one more full-domain ORION field on **every** rank — about an
eighth of the fine level per level on a 3-D mesh, a quarter on a 2-D one. The
levels are built at setup and held for the whole run, not freed as the solve
leaves them behind.

**Problem size.** Everything on this page is strong scaling of one 7 M-cell
case. Larger meshes, more families, or two-way coupling will behave differently.

---

## 12. Quick reference

| Goal | Approach |
| --- | --- |
| Cover all sockets | ≥ 1 rank per socket, or `OMP_PROC_BIND=spread` |
| Keep NUMA locality | Keep each OpenMP team inside one socket; do not interleave pages ([§6](#6-strong-scaling)) |
| Reduce memory pressure | Fewer ranks, more threads — free in speed on the pre-wave baseline; check [§6](#6-strong-scaling) before relying on it on the current tree |
| Avoid idle ranks | Cut the mesh into at least as many blocks as ranks |
| Avoid many small blocks per rank | Especially at `1 × N` — see [§11](#11-known-scalability-limits) |
| Scale to more nodes | Replicate the same per-node layout |
| Measure the solver | `shell-diter`, last window, minimum of repetitions |
| Measure efficiency | `core-cycles/iter`, not seconds |
| Validate placement | `I_MPI_DEBUG=4`, check the rank-to-core map |

---

## 13. Job index

Every job referenced on this page, its tree, its node set, and what it ran.
Results under `results/<job>.tsv`; raw per-repetition logs under
`runs/<job>/logs/`.

| Job | Tree | Nodes | What it ran |
| --- | --- | --- | --- |
| 302995 | `17bc107`+ | other | First full scaling matrix. Campaign design and harness: Marco Grossi. |
| 303350 | `aeb34a5` (`origin/merging`, before this session's fixes) + `ICE_BIN_DIR` | wn[05-07] | Baseline full matrix, 63 runs — source of every table in §6–§11 unless stated otherwise. |
| 303351 | same tree as 303350 | wn08 | OpenMP placement sweep (packed vs. spread across sockets), plus a first interleaved-pages test at `1 × 40`/`1 × 80`. |
| 303387 | `5265964` (W1, W3, W4, W5 + timer regions) | wn[08-10] | Reduced matrix, mid-wave (before C1, C9, C4). Source of the halo breakdown's "before" column. |
| 303402 | `e7cdef7` (+ C1, C9a, C9b) | wn[04-06] | Reduced matrix. Source of the halo breakdown's "after" column. |
| 303417 | `5575f65` (+ C4a, C4b) | wn[05-07], same nodes as 303350 | Reduced matrix: the wave's before/after at the four points in §6, plus `blk-C80-1x80-R24`. |
| 303418 | `5575f65`, same binary as 303417 | wn09 | NUMA: plain (first-touch) vs. interleaved pages at `1 × 40` and `1 × 80`. |
| 303439 | `3d43d0e` (+ C5b, slab tiles) | wn[05-07] | Reduced matrix: C5b neutral at 1 × 80 (0.338 s), −31/−40 % on the 20-thread socket legs, many-blocks leg doubled → C5b reverted (`461eb4f`). |
| 303482 | `461eb4f` (+ C3, parallel first touch) | wn[05-07] | Reduced matrix: C3 changed nothing (1 × 80 0.3373 s) → reverted. |
| 303484 | `461eb4f` (the final solver code plus C3, measured neutral and reverted afterwards) | 3 nodes | Greedy vs. halo-objective cuts (ATLAS MDB) at R12, R16, R20, R24 — §6. (303478, its first submission, ran the reverted C5b binary and was cancelled.) |
| 303485 | `461eb4f` (the final solver code plus C3, measured neutral and reverted afterwards) | wn[05-07], same nodes as 303350 | The full nine-point matrix on the final tree, 63 runs — the second table of §6 and the before/after of every leg. (303481, its first submission, was cancelled for the same reason as 303478.) |
