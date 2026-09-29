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
feed threads is never necessary.

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
| 80 | `1 × 80` | fastest measured, lowest memory, needs explicit placement |
| 160 | `8 × 20` | two nodes, same per-node layout |
| 240 | `12 × 20` | three nodes, same per-node layout |

### Configurations to avoid

* An OpenMP team packed onto **one socket** while other sockets sit idle.
* More ranks than blocks.
* A rank count whose memory exceeds the node.

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
the reference benchmark it prints:

```text
 ICE Timing | Iter 40 | 10 iters | wall/iter  3.2827E-01 s | rank min ... avg ...
            | imbalance  0.0 % | exchange wait  7.5 % | collective wait  2.9 %
 ICE Phases | source  3.8175E-02 s ( 11.5 %) | flux  2.3502E-01 s ( 70.7 %)
            | halo  4.9976E-02 s ( 15.0 %)
 ICE Ranks  | compute/iter max ... | min ... | mean ... | spread   5.7 %
 ICE Cycles | core-cycles/iter  6.8337E+10 | aggregate    205.60 GHz
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
to 3.83 GHz while eighty sit at 2.57. That factor of 0.67 is the chip's power
management, not ICE, and it accounts for essentially the whole difference
between ICE's wall-clock efficiency at a full node (54 %) and its core-cycle
efficiency (80 %).

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
physical cores per node and no hyperthreading.

<figure>
  {% include "user/images/scaling-envelope.svg" %}
</figure>

| Cores | Nodes | MPI × OpenMP | Wall/iter | Efficiency | Ghost-credited |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | `1 × 1` | 14.262 s | 100.0 % | 100.0 % |
| 4 | 1 | `4 × 1` | 3.640 s | 97.2 % | 101.3 % |
| 8 | 1 | `4 × 2` | 1.887 s | 94.6 % | 98.6 % |
| 16 | 1 | `4 × 4` | 0.988 s | 94.6 % | 98.6 % |
| 20 | 1 | `4 × 5` | 0.809 s | 93.7 % | 97.7 % |
| 40 | 1 | `4 × 10` | 0.472 s | 87.8 % | 91.5 % |
| 80 | 1 | `4 × 20` | 0.328 s | 80.0 % | 83.4 % |
| 160 | 2 | `8 × 20` | 0.164 s | 79.1 % | 84.0 % |
| 240 | 3 | `12 × 20` | 0.114 s | 75.3 % | 83.1 % |

Efficiency stays above **93 % through a quarter node**, reaches **80 % at a full
node**, and is then flat across two and three nodes.

**The two efficiency columns differ because the decomposition is real work.**
Each rank count runs its own cut of the mesh, and cutting adds ghost cells:
4.2 % more cells at 4 ranks, 10.4 % at 12. The raw column charges ICE for them;
the ghost-credited column does not. Quote the raw column when asking "what does
this many ranks cost me", and the ghost-credited column when asking "how well
does the parallelisation itself hold up".

---

## 7. MPI, OpenMP and hybrid execution

<figure>
  {% include "user/images/scaling-modes.svg" %}
</figure>

| Cores | Pure OpenMP `1 × N` | Pure MPI `N × 1` | Hybrid `4 × N/4` |
| ---: | ---: | ---: | ---: |
| 4 | **99.0 %** | 97.2 % | 97.2 % |
| 8 | **96.4 %** | 93.9 % | 94.6 % |
| 16 | **97.1 %** | 89.7 % | 94.6 % |
| 20 | **97.1 %** | — | 93.7 % |
| 24 | — | 68.5 % | — |
| 40 | **91.4 %** | — | 87.8 % |
| 80 | **83.4 %** | — | 80.0 % |

All three models work, and the spread between them is small compared with what
placement does. Pure OpenMP reads highest because it runs an undecomposed mesh
and therefore computes no ghost cells at all; per cell of real physics the
hybrid and MPI runs are the more efficient.

Three things to read carefully:

1. **The pure modes are pinned across all four sockets at every core count**, so
   they are compared on equal hardware rather than on whatever the default
   placement gives. The hybrid points are launched the way a user would, with
   `I_MPI_PIN_DOMAIN=omp` and no mask, and below a full node those domains are
   assigned in core order and need not reach every socket. At 20 cores the same
   `4 × 5` layout with an explicit spread mask reads 95.1 % against the table's
   93.7 %, which bounds the effect at about **1.4 points**. Only the 80-core
   points fill the node and are free of it entirely.
2. **Pure OpenMP's advantage is the ghost cells, not the threading.** At 80
   cores it computes 4.2–12.5 % fewer cells than any of the hybrid layouts.
   Credit those cells back and the ranking inverts — see
   [§8](#8-choosing-the-rankthread-split).
3. **Pure MPI's 24-rank point carries 18.8 % more cells** than the reference.
   Credited for them it reads 81.3 %, in line with the rest. What ends the
   pure-MPI curve is memory ([§10](#10-memory)), not efficiency.

---

## 8. Choosing the rank/thread split

<figure>
  {% include "user/images/scaling-layout.svg" %}
</figure>

At 80 cores on one node:

| Layout | Ranks per socket | Wall/iter | vs best | Node memory |
| --- | ---: | ---: | ---: | ---: |
| `1 × 80` | — | 0.3184 s | best | 5.7 GB |
| `8 × 10` | 2 | 0.3238 s | +1.7 % | 35.6 GB |
| `4 × 20` | 1 | 0.3284 s | +3.1 % | 20.1 GB |
| `16 × 5` | 4 | 0.3334 s | +4.7 % | 72.7 GB |

**Do not over-tune this.** The whole family spans 4.7 %, which is close to the
measurement floor, while the memory it costs spans a factor of thirteen. Choose
the rank count on memory and on the block count, and move on.

---

## 9. Scaling across nodes

<figure>
  {% include "user/images/scaling-nodes.svg" %}
</figure>

Adding nodes is nearly free. Referred to the 80-core point on one node:

| | Speed-up | Of linear |
| --- | ---: | ---: |
| 160 cores, 2 nodes | 1.998 × | **99.9 %** |
| 240 cores, 3 nodes | 2.881 × | 96.0 % |

The first node boundary costs nothing measurable. Splitting the iteration by
phase locates the remaining 4 %: from 80 to 240 cores the **compute scales at
103 % of linear and the halo at 69 %**. The limit is ICE's own ghost fill and
boundary work, which does not get cheaper as the messages get smaller — not the
interconnect. `exchange wait` is identical at one node and two, and only moves
at three.

The practical consequence: **if wall-clock time is what matters, spreading a
fixed case over more nodes is close to free in efficiency.** It costs node-hours
in proportion, so it is a trade against allocation budget, not against
scalability.

---

## 10. Memory

<figure>
  {% include "user/images/scaling-memory.svg" %}
</figure>

Every MPI rank allocates the **whole domain** and updates only its own blocks,
so the node footprint grows linearly with rank count instead of shrinking with
the decomposition:

| Ranks | 1 | 2 | 4 | 8 | 16 | 24 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GB/node | 5.7 | 10.7 | 20.1 | 35.6 | 72.7 | 112.4 |

That is **~4.7 GB per rank** at 7 M cells, whatever the rank count. A 187 GB
node therefore holds about **40 ranks** at this problem size, and pure MPI at a
full node is impossible: `80 × 1` would need some 375 GB.

This is why the rank count should be chosen on memory. It is also why more
OpenMP threads per rank is the right answer when memory is tight — and, on the
current code, it costs nothing in speed to do so.

---

## 11. Known scalability limits

**Memory per rank.** Domain-wide structures live on every rank. Very large rank
counts exhaust memory before they exhaust cores.

**Block count.** Ranks beyond the number of blocks idle, and an uneven split
shows up directly in the `balance` figure at startup. Cut the mesh for the rank
count you intend to use.

**Cell balance is not work balance.** At reported 100 % cell-count balance
at every rank count in this benchmark, ICE still measures up to 8 % spread
in per-rank compute time. Two reasons: the reconstruction retries on unphysical
states, so per-cell cost is data-dependent; and a rank's mix of physical
boundaries and cut connections is not balanced by balancing cells — a symmetry
condition costs far more than a connection's straight copy.

**The halo phase.** At 240 cores it is 21 % of the iteration and scales at 69 %
of linear. It is the limit on multi-node runs.

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
| Keep NUMA locality | Keep each OpenMP team inside one socket |
| Reduce memory pressure | Fewer ranks, more threads — it costs nothing in speed |
| Avoid idle ranks | Cut the mesh into at least as many blocks as ranks |
| Scale to more nodes | Replicate the same per-node layout |
| Measure the solver | `shell-diter`, last window, minimum of repetitions |
| Measure efficiency | `core-cycles/iter`, not seconds |
| Validate placement | `I_MPI_DEBUG=4`, check the rank-to-core map |

