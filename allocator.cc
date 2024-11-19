#include <linux/mman.h>
#include <sys/mman.h>
#include <unistd.h>
#include <map>
#include <iostream>
#include <cstring>
#include <numa.h>

#include "allocator.h"
#include "macros.h"
#include "sh_macros.h"
#include "spinlock.h"
#include "lockguard.h"
#include "static_vector.h"
#include "counter.h"
#ifdef SHMM_CACHE_OPT
extern "C" {
#include "ulcc/ulcc.h"
}
#endif

using namespace util;

static event_counter evt_allocator_total_region_usage(
    "allocator_total_region_usage_bytes");

// page+alloc routines taken from masstree

#ifdef MEMCHECK_MAGIC
const allocator::pgmetadata *
allocator::PointerToPgMetadata(const void *p)
{
  static const size_t hugepgsize = GetHugepageSize();
  if (unlikely(!ManagesPointer(p)))
    return nullptr;
  const size_t cpu = PointerToCpu(p);
  const regionctx &pc = g_regions[cpu];
  if (p >= pc.region_begin)
    return nullptr;
  // round pg down to page
  p = (const void *) ((uintptr_t)p & ~(hugepgsize-1));
  const pgmetadata *pmd = (const pgmetadata *) p;
  ALWAYS_ASSERT((pmd->unit_ % AllocAlignment) == 0);
  ALWAYS_ASSERT((MAX_ARENAS * AllocAlignment) >= pmd->unit_);
  return pmd;
}
#endif

size_t
allocator::GetHugepageSizeImpl()
{
  FILE *f = fopen("/proc/meminfo", "r");
  assert(f);
  char *linep = NULL;
  size_t n = 0;
  static const char *key = "Hugepagesize:";
  static const int keylen = strlen(key);
  size_t size = 0;
  while (getline(&linep, &n, f) > 0) {
    if (strstr(linep, key) != linep)
      continue;
    size = atol(linep + keylen) * 1024;
    break;
  }
  fclose(f);
  assert(size);
  return size;
}

size_t
allocator::GetPageSizeImpl()
{
  return sysconf(_SC_PAGESIZE);
}

bool
allocator::UseMAdvWillNeed()
{
  static const char *px = getenv("DISABLE_MADV_WILLNEED");
  static const std::string s = px ? to_lower(px) : "";
  static const bool use_madv = !(s == "1" || s == "true");
  return use_madv;
}

void
allocator::Initialize(size_t ncpus, size_t maxpercore)
{
  static spinlock s_lock;
  static bool s_init = false;
  if (likely(s_init))
    return;
  lock_guard<spinlock> l(s_lock);
  if (s_init)
    return;
  ALWAYS_ASSERT(!g_memstart);
  ALWAYS_ASSERT(!g_memend);
  ALWAYS_ASSERT(!g_ncpus);
  ALWAYS_ASSERT(!g_maxpercore);

  static const size_t hugepgsize = GetHugepageSize();

  // round maxpercore to the nearest hugepagesize
  maxpercore = slow_round_up(maxpercore, hugepgsize);

  g_ncpus = ncpus;
  g_maxpercore = maxpercore;

  // mmap() the entire region for now, but just as a marker
  // (this does not actually cause physical pages to be allocated)
  // note: we allocate an extra hugepgsize so we can guarantee alignment
  // of g_memstart to a huge page boundary

  void * const x = mmap(nullptr, g_ncpus * g_maxpercore + hugepgsize,
      PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (x == MAP_FAILED) {
    perror("mmap");
    ALWAYS_ASSERT(false);
  }

  void * const endpx = (void *) ((uintptr_t)x + g_ncpus * g_maxpercore + hugepgsize);
  std::cerr << "allocator::Initialize()" << std::endl
            << "  hugepgsize: " << hugepgsize << std::endl
            << "  use MADV_WILLNEED: " << UseMAdvWillNeed() << std::endl
            << "  mmap() region [" << x << ", " << endpx << ")" << std::endl;

  g_memstart = reinterpret_cast<void *>(util::iceil(uintptr_t(x), hugepgsize));
  g_memend = reinterpret_cast<char *>(g_memstart) + (g_ncpus * g_maxpercore);

  ALWAYS_ASSERT(!(reinterpret_cast<uintptr_t>(g_memstart) % hugepgsize));
  ALWAYS_ASSERT(reinterpret_cast<uintptr_t>(g_memend) <=
      (reinterpret_cast<uintptr_t>(x) + (g_ncpus * g_maxpercore + hugepgsize)));

  for (size_t i = 0; i < g_ncpus; i++) {
    g_regions[i].region_begin =
      reinterpret_cast<char *>(g_memstart) + (i * g_maxpercore);
    g_regions[i].region_end   =
      reinterpret_cast<char *>(g_memstart) + ((i + 1) * g_maxpercore);
    std::cerr << "cpu" << i << " owns [" << g_regions[i].region_begin
              << ", " << g_regions[i].region_end << ")" << std::endl;
    ALWAYS_ASSERT(g_regions[i].region_begin < g_regions[i].region_end);
    ALWAYS_ASSERT(g_regions[i].region_begin >= x);
    ALWAYS_ASSERT(g_regions[i].region_end <= endpx);
  }

  s_init = true;
}

void allocator::InitializeWorkerMemory(size_t nworkers)
{
  ALWAYS_ASSERT(!g_worker_memstart);
  ALWAYS_ASSERT(!g_worker_memend);
  ALWAYS_ASSERT(!g_nworkers);
  ALWAYS_ASSERT(!g_maxperworker);

  static const size_t pgsize = 4096;
  g_maxperworker = sizeof(worker_numa_memory_layout) + g_mem_per_worker;
  g_maxperworker = slow_round_up(g_maxperworker, pgsize);
  g_nworkers = nworkers;

  void *const x = mmap(nullptr, g_nworkers * g_maxperworker,
          PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (x == MAP_FAILED) {
    perror("mmap");
    ALWAYS_ASSERT(false);
  }
  void *const endpx = (void *) ((uintptr_t)x + g_nworkers * g_maxperworker);
  std::cerr << "allocator::InitializeWorkerMemory()" << std::endl
            << "  pgsize: " << pgsize << std::endl
            //<< "  use MADV_WILLNEED: " << UseMAdvWillNeed() << std::endl
            << "  mmap() region [" << x << ", " << endpx << ")" << std::endl;

  g_worker_memstart = x; //reinterpret_cast<void *>(util::iceil(uintptr_t(x), pgsize));
  g_worker_memend = endpx; //reinterpret_cast<char *>(g_worker_memstart) + (g_nworkers * g_maxperworker);
#ifdef USE_BCC
  ComputeSockConfiguration();
#endif
}

void
allocator::DumpStats()
{
  std::cerr << "[allocator] ncpus=" << g_ncpus << std::endl;
  for (size_t i = 0; i < g_ncpus; i++) {
    const bool f = g_regions[i].region_faulted;
    const size_t remaining =
      intptr_t(g_regions[i].region_end) -
      intptr_t(g_regions[i].region_begin);
    std::cerr << "[allocator] cpu=" << i << " fully_faulted?=" << f
              << " remaining=" << remaining << " bytes" << std::endl;
  }
}

static void *
initialize_page(void *page, const size_t pagesize, const size_t unit)
{
  INVARIANT(((uintptr_t)page % pagesize) == 0);

#ifdef MEMCHECK_MAGIC
  ::allocator::pgmetadata *pmd = (::allocator::pgmetadata *) page;
  pmd->unit_ = unit;
  page = (void *) ((uintptr_t)page + sizeof(*pmd));
#endif

  void *first = (void *)util::iceil((uintptr_t)page, (uintptr_t)unit);
  INVARIANT((uintptr_t)first + unit <= (uintptr_t)page + pagesize);
  void **p = (void **)first;
  void *next = (void *)((uintptr_t)p + unit);
  while ((uintptr_t)next + unit <= (uintptr_t)page + pagesize) {
    INVARIANT(((uintptr_t)p % unit) == 0);
    *p = next;
#ifdef MEMCHECK_MAGIC
    NDB_MEMSET(
        (char *) p + sizeof(void **),
        MEMCHECK_MAGIC, unit - sizeof(void **));
#endif
    p = (void **)next;
    next = (void *)((uintptr_t)next + unit);
  }
  INVARIANT(((uintptr_t)p % unit) == 0);
  *p = NULL;
#ifdef MEMCHECK_MAGIC
  NDB_MEMSET(
      (char *) p + sizeof(void **),
      MEMCHECK_MAGIC, unit - sizeof(void **));
#endif
  return first;
}

void *
allocator::AllocateArenas(size_t cpu, size_t arena)
{
  INVARIANT(cpu < g_ncpus);
  INVARIANT(arena < MAX_ARENAS);
  INVARIANT(g_memstart);
  INVARIANT(g_maxpercore);
  static const size_t hugepgsize = GetHugepageSize();

  regionctx &pc = g_regions[cpu];
  pc.lock.lock();
  if (likely(pc.arenas[arena])) {
    // claim
    void *ret = pc.arenas[arena];
    pc.arenas[arena] = nullptr;
    pc.lock.unlock();
    return ret;
  }

  void * const mypx = AllocateUnmanagedWithLock(pc, 1); // releases lock
  return initialize_page(mypx, hugepgsize, (arena + 1) * AllocAlignment);
}

void *
allocator::AllocateUnmanaged(size_t cpu, size_t nhugepgs)
{
  regionctx &pc = g_regions[cpu];
  pc.lock.lock();
  return AllocateUnmanagedWithLock(pc, nhugepgs); // releases lock
}

void *
allocator::AllocateUnmanagedWithLock(regionctx &pc, size_t nhugepgs)
{
  static const size_t hugepgsize = GetHugepageSize();

  void * const mypx = pc.region_begin;

  // check alignment
  if (reinterpret_cast<uintptr_t>(mypx) % hugepgsize)
    ALWAYS_ASSERT(false);

  void * const mynewpx =
    reinterpret_cast<char *>(mypx) + nhugepgs * hugepgsize;

  if (unlikely(mynewpx > pc.region_end)) {
    std::cerr << "allocator::AllocateUnmanagedWithLock():" << std::endl
              << "  region ending at " << pc.region_end << " OOM" << std::endl;
    ALWAYS_ASSERT(false); // out of memory otherwise
  }

  const bool needs_mmap = !pc.region_faulted;
  pc.region_begin = mynewpx;
  pc.lock.unlock();

  evt_allocator_total_region_usage.inc(nhugepgs * hugepgsize);

  if (needs_mmap) {
    void * const x = mmap(mypx, hugepgsize, PROT_READ | PROT_WRITE,
        MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED, -1, 0);
    if (unlikely(x == MAP_FAILED)) {
      perror("mmap");
      ALWAYS_ASSERT(false);
    }
    INVARIANT(x == mypx);
    const int advice =
      UseMAdvWillNeed() ? MADV_HUGEPAGE | MADV_WILLNEED : MADV_HUGEPAGE;
    if (madvise(x, hugepgsize, advice)) {
      perror("madvise");
      ALWAYS_ASSERT(false);
    }
  }

  return mypx;
}

void
allocator::ReleaseArenas(void **arenas)
{
  // cpu -> [(head, tail)]
  // XXX: use a small_map here?
  std::map<size_t, static_vector<std::pair<void *, void *>, MAX_ARENAS>> m;
  for (size_t arena = 0; arena < MAX_ARENAS; arena++) {
    void *p = arenas[arena];
    while (p) {
      void * const pnext = *reinterpret_cast<void **>(p);
      const size_t cpu = PointerToCpu(p);
      auto it = m.find(cpu);
      if (it == m.end()) {
        auto &v = m[cpu];
        v.resize(MAX_ARENAS);
        *reinterpret_cast<void **>(p) = nullptr;
        v[arena].first = v[arena].second = p;
      } else {
        auto &v = it->second;
        if (!v[arena].second) {
          *reinterpret_cast<void **>(p) = nullptr;
          v[arena].first = v[arena].second = p;
        } else {
          *reinterpret_cast<void **>(p) = v[arena].first;
          v[arena].first = p;
        }
      }
      p = pnext;
    }
  }
  for (auto &p : m) {
    INVARIANT(!p.second.empty());
    regionctx &pc = g_regions[p.first];
    lock_guard<spinlock> l(pc.lock);
    for (size_t arena = 0; arena < MAX_ARENAS; arena++) {
      INVARIANT(bool(p.second[arena].first) == bool(p.second[arena].second));
      if (!p.second[arena].first)
        continue;
      *reinterpret_cast<void **>(p.second[arena].second) = pc.arenas[arena];
      pc.arenas[arena] = p.second[arena].first;
    }
  }
}

static void
numa_hint_memory_placement(void *px, size_t sz, unsigned node)
{
  struct bitmask *bm = numa_allocate_nodemask();
  numa_bitmask_setbit(bm, node);
  numa_interleave_memory(px, sz, bm);
  numa_free_nodemask(bm);
}

void
allocator::FaultRegion(size_t cpu)
{
  static const size_t hugepgsize = GetHugepageSize();
  ALWAYS_ASSERT(cpu < g_ncpus);
  regionctx &pc = g_regions[cpu];
  if (pc.region_faulted)
    return;
  lock_guard<std::mutex> l1(pc.fault_lock);
  lock_guard<spinlock> l(pc.lock); // exclude other users of the allocator
  if (pc.region_faulted)
    return;
  // mmap the entire region + memset it for faulting
  if (reinterpret_cast<uintptr_t>(pc.region_begin) % hugepgsize)
    ALWAYS_ASSERT(false);
  const size_t sz =
    reinterpret_cast<uintptr_t>(pc.region_end) -
    reinterpret_cast<uintptr_t>(pc.region_begin);
  void * const x = mmap(pc.region_begin, sz, PROT_READ | PROT_WRITE,
      MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED, -1, 0);
  if (unlikely(x == MAP_FAILED)) {
    perror("mmap");
    std::cerr << "  cpu" << cpu
              << " [" << pc.region_begin << ", " << pc.region_end << ")"
              << std::endl;
    ALWAYS_ASSERT(false);
  }
  ALWAYS_ASSERT(x == pc.region_begin);
  const int advice =
    UseMAdvWillNeed() ? MADV_HUGEPAGE | MADV_WILLNEED : MADV_HUGEPAGE;
  if (madvise(x, sz, advice)) {
    perror("madvise");
    ALWAYS_ASSERT(false);
  }
  numa_hint_memory_placement(
      pc.region_begin,
      (uintptr_t)pc.region_end - (uintptr_t)pc.region_begin,
      numa_node_of_cpu(cpu));
  const size_t nfaults =
    ((uintptr_t)pc.region_end - (uintptr_t)pc.region_begin) / hugepgsize;
  std::cerr << "cpu" << cpu << " starting faulting region ("
            << intptr_t(pc.region_end) - intptr_t(pc.region_begin)
            << " bytes / " << nfaults << " hugepgs)" << std::endl;
  timer t;
  for (char *px = (char *) pc.region_begin;
       px < (char *) pc.region_end;
       px += CACHELINE_SIZE)
    *px = 0xDE;
  std::cerr << "cpu" << cpu << " finished faulting region in "
            << t.lap_ms() << " ms" << std::endl;
  pc.region_faulted = true;
}

void allocator::FaultWorkerMemory(size_t wid, size_t pin_cpu)
{
  ALWAYS_ASSERT(wid < g_nworkers);
  //static const size_t hugepgsize = GetHugepageSize();
  void *mem_start = reinterpret_cast<char *>(g_worker_memstart) +
          (wid * g_maxperworker);
  void *mem_end = reinterpret_cast<char *>(mem_start) + g_maxperworker;

  // mmap the entire region + memset it for faulting
  void *const x = mmap(mem_start, g_maxperworker, PROT_READ | PROT_WRITE,
      MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED, -1, 0);
  if (unlikely(x == MAP_FAILED)) {
    perror("mmap");
    std::cerr << "  worker" << wid
              << " [" << mem_start << ", " << mem_end << ")"
              << std::endl;
    ALWAYS_ASSERT(false);
  }
  ALWAYS_ASSERT(x == mem_start);
//  const int advice =
//    UseMAdvWillNeed() ? MADV_HUGEPAGE | MADV_WILLNEED : MADV_HUGEPAGE;
//  if (madvise(x, g_maxperworker, advice)) {
//    perror("madvise");
//    ALWAYS_ASSERT(false);
//  }
  numa_hint_memory_placement(mem_start, g_maxperworker,
      numa_node_of_cpu(pin_cpu));
  //const size_t nfaults = g_maxperworker / hugepgsize;
  //std::cerr << "worker" << wid << " starting faulting worker memory ("
  //          << g_maxperworker
  //          << " bytes / " << nfaults << " hugepgs)" << std::endl;
  //timer t;
  for (char *px = reinterpret_cast<char *>(mem_start);
          px < reinterpret_cast<char *>(mem_end); px += 4096)
    *px = 'x';

#ifdef SHMM_CACHE_OPT
  const int strong_colors_per_worker =
      static_cast<int>(LLC_COLORS_STRONG_LOCALITY / g_nworkers);
  ALWAYS_ASSERT(strong_colors_per_worker > 0);
  const int strong_color_low = static_cast<int>(strong_colors_per_worker * wid);
  const int strong_color_high = strong_color_low + strong_colors_per_worker - 1;
  const int weak_color_low = LLC_COLORS_STRONG_LOCALITY;
  const int weak_color_high = weak_color_low + LLC_COLORS_WEAK_LOCALITY - 1;
  // Remapping mempool 0 to strong locality region
  cc_cacheregn_t regn;
  cc_cacheregn_clr(&regn);
  cc_cacheregn_set(&regn, strong_color_low, strong_color_high, 1);
  unsigned long remap_start = reinterpret_cast<unsigned long &>(mem_start);
  remap_start += sizeof(worker_numa_memory_layout);
  unsigned long remap_end = remap_start + SHMM_MEMPOOL_SIZE_0;
  int remap_ret = cc_remap(&remap_start, &remap_end, 1, &regn, 0);
  ALWAYS_ASSERT(remap_ret != -1);
  // Remapping mempool 1 to weak locality region
  cc_cacheregn_clr(&regn);
  cc_cacheregn_set(&regn, weak_color_low, weak_color_high, 1);
  remap_start = remap_end;
  remap_end += SHMM_MEMPOOL_SIZE_1;
  remap_ret = cc_remap(&remap_start, &remap_end, 1, &regn, 0);
  ALWAYS_ASSERT(remap_ret != -1);
#endif
  //std::cerr << "worker" << wid << " finished faulting worker memory in "
  //          << t.lap_ms() << " ms" << std::endl;
}

#ifdef USE_BCC
void allocator::ComputeSockConfiguration()
{
  ComputeSockConfiguration2();
  //PrintSockConfiguration();
}

void allocator::PrintSockConfiguration()
{
  std::cerr << "nsocks_used = " << nsocks_used << std::endl;

  std::cerr << "nworkers_on_sock:";
  for (size_t sockid = 0; sockid < nsocks_used; ++sockid)
    std::cerr << " " << nworkers_on_sock[sockid];
  std::cerr << std::endl;

  std::cerr << "wid_to_coreid_mapping:";
  for (size_t wid = 0; wid < g_nworkers; ++wid)
    std::cerr << " " << wid_to_coreid_mapping[wid];
  std::cerr << std::endl;

  std::cerr << "wid_to_sockid_mapping:";
  for (size_t wid = 0; wid < g_nworkers; ++wid)
    std::cerr << " " << wid_to_sockid_mapping[wid];
  std::cerr << std::endl;

  std::cerr << "wid_to_idx_in_sock_mapping:";
  for (size_t wid = 0; wid < g_nworkers; ++wid)
    std::cerr << " " << wid_to_idx_in_sock_mapping[wid];
  std::cerr << std::endl;

  std::cerr << "wid_to_idx_in_tids_read:";
  for (size_t wid = 0; wid < g_nworkers; ++wid)
    std::cerr << " " << wid_to_idx_in_tids_read[wid];
  std::cerr << std::endl;
}

void allocator::ComputeSockConfiguration1()
{
  ALWAYS_ASSERT(g_nworkers > 0);
  nsocks_used = ((g_nworkers + 7) & (~7)) >> 3;
  size_t nworkers_per_sock = 0, nworkers_last_sock = 0;
  if (nsocks_used == 1) {
    nworkers_per_sock = g_nworkers;
    nworkers_last_sock = g_nworkers;
  } else if ((g_nworkers % 8) == 0) {
    nworkers_per_sock = 8;
    nworkers_last_sock = 8;
  } else {
    nworkers_per_sock = 8;
    nworkers_last_sock = g_nworkers % 8;
    while (nworkers_per_sock - 1 >= nworkers_last_sock + nsocks_used - 1) {
      --nworkers_per_sock;
      nworkers_last_sock += nsocks_used - 1;
    }
  }
  ALWAYS_ASSERT(nsocks_used > 0 && nworkers_per_sock > 0 &&
      nworkers_last_sock > 0);
  for (size_t isock = 0; isock < nsocks_used - 1; ++isock)
    nworkers_on_sock[isock] = nworkers_per_sock;
  nworkers_on_sock[nsocks_used - 1] = nworkers_last_sock;
  for (size_t wid = 0; wid < g_nworkers; ++wid) {
    size_t sockid = wid / nworkers_per_sock;
    size_t idx_in_sock = wid % nworkers_per_sock;
    wid_to_sockid_mapping[wid] = sockid;
    wid_to_idx_in_sock_mapping[wid] = idx_in_sock;
    wid_to_idx_in_tids_read[wid] = wid;
    wid_to_coreid_mapping[wid] = coreids_mapping[sockid][idx_in_sock];
  }
}

void allocator::ComputeSockConfiguration2()
{
  ALWAYS_ASSERT(g_nworkers > 0);
  nsocks_used = ((g_nworkers + 7) & (~7)) >> 3;
  //nsocks_used = (g_nworkers <= 32) ? 4 : 8;
  for (size_t isock = 0; isock < g_nworkers % nsocks_used; ++isock)
    nworkers_on_sock[isock] = g_nworkers / nsocks_used + 1;
  for (size_t isock = g_nworkers % nsocks_used; isock < nsocks_used; ++isock)
    nworkers_on_sock[isock] = g_nworkers / nsocks_used;
  for (size_t wid = 0; wid < g_nworkers; ++wid) {
    size_t sockid = wid % nsocks_used;
    size_t idx_in_sock = wid / nsocks_used;
    wid_to_sockid_mapping[wid] = sockid;
    wid_to_idx_in_sock_mapping[wid] = idx_in_sock;
    wid_to_coreid_mapping[wid] = coreids_mapping[sockid][idx_in_sock];
  }
  size_t idx_in_tids_read = 0;
  for (size_t isock = 0; isock < nsocks_used; ++isock) {
    for (size_t iw = 0; iw < nworkers_on_sock[isock]; ++iw) {
      size_t wid = iw * nsocks_used + isock;
      ALWAYS_ASSERT(wid < g_nworkers);
      wid_to_idx_in_tids_read[wid] = idx_in_tids_read++;
    }
  }
}
#endif

void *allocator::g_memstart = nullptr;
void *allocator::g_memend = nullptr;
size_t allocator::g_ncpus = 0;
size_t allocator::g_maxpercore = 0;
percore<allocator::regionctx> allocator::g_regions;

void *allocator::g_worker_memstart = nullptr;
void *allocator::g_worker_memend = nullptr;
size_t allocator::g_nworkers = 0;
size_t allocator::g_maxperworker = 0;
const size_t allocator::g_mem_per_worker = BCC_PER_WORKER_MEMORY;

//XXX: The mapping needs to be adjusted for different platforms.
const size_t allocator::coreids_mapping[8][8] = {
    {0, 33, 34, 35, 36, 37, 38, 39},
    {1, 2, 3, 4, 5, 6, 7, 8},
    {9, 10, 11, 12, 13, 14, 15, 16},
    {17, 18, 19, 20, 21, 22, 23, 24},
    {25, 26, 27, 28, 29, 30, 31, 32},
    {40, 41, 42, 43, 44, 45, 46, 47},
    {48, 49, 50, 51, 52, 53, 54, 55},
    {56, 57, 58, 59, 60, 61, 62, 63}
};

#ifdef USE_BCC
size_t allocator::wid_to_coreid_mapping[BCC_MAX_NTHREADS] = {0,};
size_t allocator::wid_to_sockid_mapping[BCC_MAX_NTHREADS] = {0,};
size_t allocator::wid_to_idx_in_sock_mapping[BCC_MAX_NTHREADS] = {0,};
size_t allocator::wid_to_idx_in_tids_read[BCC_MAX_NTHREADS] = {0,};
size_t allocator::nworkers_on_sock[8] = {0,};
size_t allocator::nsocks_used = 0;
#endif
