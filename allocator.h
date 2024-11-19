#ifndef _NDB_ALLOCATOR_H_
#define _NDB_ALLOCATOR_H_

#include <cstdint>
#include <iterator>
#include <mutex>

#include "sh_macros.h"
#include "sh_mm.h"
#include "sh_txn.h"
#include "util.h"
#include "core.h"
#include "macros.h"
#include "spinlock.h"

// The layout of numa memory allocated for each worker thread.
struct worker_numa_memory_layout {
    // tids of all workers when a txn starts up.
    uint64_t tids_on_txn_startup[BCC_MAX_NTHREADS];
    // tids of all workers when a txn enters validation.
    uint64_t tids_on_txn_validation[BCC_MAX_NTHREADS];
    // Whether tids_on_txn_validation contains valid tids data from
    // previous txn execution.
    bool prev_tids_on_txn_validation_valid;
    // Pointer to the TID array on each socket. tids_on_txn_startup
    // and tids_on_txn_validation will be populated by copying tids
    // from tid arrays pointed to by tid_arrays.
    uint64_t *tid_arrays[8];
    // Temporary buffer for storing unique non-insert tuples in the
    // write set of each txn.
    const void *unique_non_inserts[BCC_MAX_UNIQUE_NON_INSERTS] __attribute__((aligned(16)));
    // TID array on the local socket of this worker thread. Only one
    // local_tid_array (that of the worker with smallest wid on a
    // socket) on each socket will be used and shared by all workers
    // on the socket.
    // NOTE: local_tid_array is read and modified by all workers on
    // a socket, and thus should reside on a separate cache line.
    uint64_t local_tid_array[8] CACHE_ALIGNED;
    // The shmm structure for managing local memory.
    shmm mm CACHE_ALIGNED;
    // The txn structure for managing committed/inflight transactions.
    shtxn txn;
    // Pointers to all workers' txn structures.
    shtxn *txns_ar[BCC_MAX_NTHREADS];
    // The start of local memory used by hash tables.
    char mem[0] __attribute__((aligned(4096))); //CACHE_ALIGNED;
};

class allocator {
public:

  // our allocator doesn't let allocations exceed maxpercore over a single core
  //
  // Initialize can be called many times- but only the first call has effect.
  //
  // w/o calling Initialize(), behavior for this class is undefined
  static void Initialize(size_t ncpus, size_t maxpercore);

  // Initialize numa memory for worker threads.
  static void InitializeWorkerMemory(size_t nworkers);

  static void DumpStats();

  // returns an arena linked-list
  static void *
  AllocateArenas(size_t cpu, size_t sz);

  // allocates nhugepgs * hugepagesize contiguous bytes from CPU's region and
  // returns the raw, unmanaged pointer.
  //
  // Note that memory returned from here cannot be released back to the
  // allocator, so this should only be used for data structures which live
  // throughput the duration of the system (ie log buffers)
  static void *
  AllocateUnmanaged(size_t cpu, size_t nhugepgs);

  static void
  ReleaseArenas(void **arenas);

  static const size_t LgAllocAlignment = 4; // all allocations aligned to 2^4 = 16
  static const size_t AllocAlignment = 1 << LgAllocAlignment;
  static const size_t MAX_ARENAS = 32;

  static inline std::pair<size_t, size_t>
  ArenaSize(size_t sz)
  {
    const size_t allocsz = util::round_up<size_t, LgAllocAlignment>(sz);
    const size_t arena = allocsz / AllocAlignment - 1;
    return std::make_pair(allocsz, arena);
  }

  // slow, but only needs to be called on initialization
  static void
  FaultRegion(size_t cpu);

  // Fault worker memory on local numa node.
  static void FaultWorkerMemory(size_t wid, size_t pin_cpu);

  // returns true if managed by this allocator, false otherwise
  static inline bool
  ManagesPointer(const void *p)
  {
    return p >= g_memstart && p < g_memend;
  }

  // assumes p is managed by this allocator- returns the CPU from which this pointer
  // was allocated
  static inline size_t
  PointerToCpu(const void *p)
  {
    ALWAYS_ASSERT(p >= g_memstart);
    ALWAYS_ASSERT(p < g_memend);
    const size_t ret =
      (reinterpret_cast<const char *>(p) -
       reinterpret_cast<const char *>(g_memstart)) / g_maxpercore;
    ALWAYS_ASSERT(ret < g_ncpus);
    return ret;
  }

#ifdef MEMCHECK_MAGIC
  struct pgmetadata {
    uint32_t unit_; // 0-indexed
  } PACKED;

  // returns nullptr if p is not managed, or has not been allocated yet.
  // p does not have to be properly aligned
  static const pgmetadata *
  PointerToPgMetadata(const void *p);
#endif

  static size_t
  GetPageSize()
  {
    static const size_t sz = GetPageSizeImpl();
    return sz;
  }

  static size_t
  GetHugepageSize()
  {
    static const size_t sz = GetHugepageSizeImpl();
    return sz;
  }

  static void *GetWorkerMemStart(int tindex)
  {
    return reinterpret_cast<char *>(g_worker_memstart) +
        tindex * g_maxperworker;
  }

  static size_t GetMemPerWorker() {
    return g_mem_per_worker;
  }

#ifdef USE_BCC
  static inline size_t wid_to_coreid(size_t wid) {
    return wid_to_coreid_mapping[wid];
  }
  static inline size_t wid_to_sockid(size_t wid) {
    return wid_to_sockid_mapping[wid];  // wid / nworkers_per_sock;
  }
  static inline size_t wid_to_idx_in_sock(size_t wid) {
    return wid_to_idx_in_sock_mapping[wid];
//    return wid / nsocks_used;
  }
  static inline bool is_wid_smallest_on_sock(size_t wid) {
    return wid_to_idx_in_sock(wid) == 0;
//    return wid_to_idx_in_sock(wid) == 0;
  }
#endif

private:
  static size_t GetPageSizeImpl();
  static size_t GetHugepageSizeImpl();
  static bool UseMAdvWillNeed();

  struct regionctx {
    regionctx()
      : region_begin(nullptr),
        region_end(nullptr),
        region_faulted(false)
    {
      NDB_MEMSET(arenas, 0, sizeof(arenas));
    }
    regionctx(const regionctx &) = delete;
    regionctx(regionctx &&) = delete;
    regionctx &operator=(const regionctx &) = delete;

    // set by Initialize()
    void *region_begin;
    void *region_end;

    bool region_faulted;

    spinlock lock;
    std::mutex fault_lock; // XXX: hacky
    void *arenas[MAX_ARENAS];
  };

  // assumes caller has the regionctx lock held, and
  // will release the lock.
  static void *
  AllocateUnmanagedWithLock(regionctx &pc, size_t nhugepgs);

#ifdef USE_BCC
  static void ComputeSockConfiguration();
  static void PrintSockConfiguration();
  static void ComputeSockConfiguration1();
  static void ComputeSockConfiguration2();
#endif

  // [g_memstart, g_memstart + ncpus * maxpercore) is the region of memory mmap()-ed
  static void *g_memstart;
  static void *g_memend; // g_memstart + ncpus * maxpercore
  static size_t g_ncpus;
  static size_t g_maxpercore;

  static percore<regionctx> g_regions CACHE_ALIGNED;

  // Worker numa memory.
  static void *g_worker_memstart;
  static void *g_worker_memend;
  static size_t g_nworkers;
  static size_t g_maxperworker;
  static const size_t g_mem_per_worker;

public:
  static const size_t coreids_mapping[8][8];
#ifdef USE_BCC
  static size_t wid_to_coreid_mapping[BCC_MAX_NTHREADS];
  static size_t wid_to_sockid_mapping[BCC_MAX_NTHREADS];
  static size_t wid_to_idx_in_sock_mapping[BCC_MAX_NTHREADS];
  static size_t wid_to_idx_in_tids_read[BCC_MAX_NTHREADS];
  static size_t nworkers_on_sock[8];
  static size_t nsocks_used;
#endif
};

#endif /* _NDB_ALLOCATOR_H_ */
