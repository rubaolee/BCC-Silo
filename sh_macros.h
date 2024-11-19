#ifndef __BCC_MACROS_H__
#define __BCC_MACROS_H__

#include "macros.h"

// NUMA memory per bench worker. Set based on statistics with -t16 -s16.
#define BCC_PER_WORKER_MEMORY   (8L * 1024L * 1024L)

// Number of shmm memory pools per worker
#define SHMM_MEMPOOLS   3
// The size of each shmm memory pool.
#define SHMM_MEMPOOL_SIZE_0   (20L * 1024L)
#define SHMM_MEMPOOL_SIZE_1   (400L * 1024L)
#define SHMM_MEMPOOL_SIZE_2   (BCC_PER_WORKER_MEMORY - SHMM_MEMPOOL_SIZE_0 - SHMM_MEMPOOL_SIZE_1)

// Max number of txns in txn_elem_list.
// Must be power of 2. Set based on statistics with -t16 -s16.
#define MAX_TXN_ELEMS     (16 * 1024)

// Enable LLC optimization for shmm management.
// #define SHMM_CACHE_OPT
// Number of colors for strong- and weak-locality shmm regions.
// NOTE: This setting is machine-dependent. Must be manually set.
// On mnode05, we have 256 colors for the 24576K LLC.
// The organization of core ids with respect to LLC sharing looks as follows:
// [1-8], [9-16], [17-24], [25-32], [33-39,0], [40-47], [48-55], [56-63]
#define LLC_COLORS_STRONG_LOCALITY  224
#define LLC_COLORS_WEAK_LOCALITY    32

// Max number of worker threads.
#define BCC_MAX_NTHREADS    32 

// For tuning performance only. Do not turn on unless
// knowing exactly what you are doing.
//#define BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS

// Optimizing start tid reading; use previous end_tid if feasible
//#define BCC_READ_STARTTID_OPT

// Check WR left conflicts in do_tuple_read
#define BCC_WR_LEFT_CONFLICT_CHECK_FORWARD

// Check invariants in BCC routines
#define BCC_CHECK

// Collect BCC statistics
//#define BCC_TXN_STATS
//#define BCC_SHMM_STATS
//#define BCC_SHTXN_STATS
//#define BCC_SHHT_STATS

#if defined(BCC_TXN_STATS) || (defined(USE_YY_CC) && (defined(BCC_SHTXN_STATS) || defined(BCC_SHM_STATS)))
#define BCC_STATS
#endif

// Max number of unique non-insert tutples in the write set of a txn
#define BCC_MAX_UNIQUE_NON_INSERTS    200

// Grouping unique non-insert write tuples together is definitely helpful
// to performance, but further using SSE does not seem to be much helpful.
//#define BCC_SEARCH_TUPLE_SSE

#define BCC_FIX_YCSB_NUMA

// Macros for computing time in us
#define TVAL(t)         (static_cast<double>((t).tv_sec) * 1.0e6 + static_cast<double>((t).tv_nsec) / 1.0e3)
#define TDIFF(t1, t2)   (TVAL(t2) - TVAL(t1))

#endif
