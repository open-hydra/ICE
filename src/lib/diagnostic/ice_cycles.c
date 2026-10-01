/* Unprivileged CPU-cycle counting for the ICE timers.
 *
 * Wall-clock seconds alone cannot answer "how well does this scale", because
 * the processor does not run at one frequency: a single active core turbos far
 * above the all-core clock (measured 3.88 GHz against a 2.1 GHz base on the
 * reference machine's Xeon Gold 6230). An efficiency built from seconds
 * therefore charges the solver for the chip's power management.
 *
 * Core-cycles are immune to that. Perfect parallelism keeps the TOTAL number of
 * core-cycles spent per iteration constant as cores are added, whatever
 * frequency each core happens to run at, so C(1)/C(N) is parallel efficiency
 * with the clock divided out.
 *
 * One counter is attached to each OpenMP thread, rather than one inherited
 * counter per process. perf_event_open's inherit flag only follows tasks
 * created AFTER the counter is enabled, and folds them in when they exit: the
 * OpenMP pool already exists by the time the timers start and its threads stay
 * alive to the end, so an inherited counter silently reports roughly the master
 * thread alone. Attaching per thread counts what is actually running.
 *
 * perf_event_open on the calling thread is permitted at perf_event_paranoid=2
 * provided kernel-mode counting is excluded. Every entry point degrades
 * quietly: if the counters cannot be opened the timers simply report no cycle
 * figure, and the wall-clock half of the instrumentation is unaffected.
 *
 * Ported from hydra-MF's hydra_cycles.c, where the design notes above were
 * established; see the ICE scaling campaign for how the numbers are used.
 */
#define _GNU_SOURCE
#include <string.h>
#include <stdint.h>

#if defined(__linux__)
#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <sys/ioctl.h>
#include <unistd.h>
#endif

#define ICE_CYC_MAX 1024

/* Stores fd+1, so that static zero-initialisation already means "unset".
   An explicit init pass would have to run from one of the threads that are
   concurrently registering, and its reset loop would race with their stores. */
static int fds1[ICE_CYC_MAX];
static int n_open = 0;

/* Attach a cycle counter to the CALLING thread, recorded in `slot`.
   Call once per OpenMP thread from inside a parallel region. */
void ice_cyc_open_thread(int *slot)
{
#if defined(__linux__)
    struct perf_event_attr a;
    int s = *slot;
    int fd;

    if (s < 0 || s >= ICE_CYC_MAX) return;
    if (fds1[s] != 0) return;                /* already attached */

    memset(&a, 0, sizeof(a));
    a.type           = PERF_TYPE_HARDWARE;
    a.size           = sizeof(a);
    a.config         = PERF_COUNT_HW_CPU_CYCLES;
    a.disabled       = 1;
    a.exclude_kernel = 1;   /* required at perf_event_paranoid=2 */
    a.exclude_hv     = 1;

    /* pid 0 = the calling thread, cpu -1 = wherever it is scheduled */
    fd = (int)syscall(__NR_perf_event_open, &a, 0, -1, -1, 0);
    if (fd < 0) return;

    ioctl(fd, PERF_EVENT_IOC_RESET, 0);
    ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
    fds1[s] = fd + 1;
    __sync_fetch_and_add(&n_open, 1);
#else
    (void)slot;
#endif
}

/* How many threads carry a counter. 0 means cycle counting is unavailable. */
void ice_cyc_count(int *n)
{
    *n = n_open;
}

/* Cycles since the counters were opened, summed over every attached thread. */
void ice_cyc_read_all(int64_t *total)
{
    int64_t sum = 0;
#if defined(__linux__)
    uint64_t v;
    int i;

    for (i = 0; i < ICE_CYC_MAX; i++) {
        if (fds1[i] == 0) continue;
        v = 0;
        if (read(fds1[i] - 1, &v, sizeof(v)) == (ssize_t)sizeof(v))
            sum += (int64_t)v;
    }
#endif
    *total = sum;
}
