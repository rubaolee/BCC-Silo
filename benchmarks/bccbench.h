#ifndef _NDB_BENCH_H_
#define _NDB_BENCH_H_

#include <climits>
#include <cstring>
#include <numa.h>
#include <stdint.h>

#include <map>
#include <numa.h>
#include <vector>
#include <utility>
#include <string>

#include "../allocator.h"
#include "../macros.h"
#include "../thread.h"
#include "../util.h"
#include "../spinbarrier.h"
#include "../rcu.h"
#include "../sh_macros.h"
#include "../sh_mm.h"
#include "../sh_txn.h"

#include "bccndb_wrapper.h"

extern void bcctest(ndb_wrapper *db, size_t nthreads, int, char**);

enum {
  RUNMODE_TIME = 0,
  RUNMODE_OPS  = 1
};

// benchmark global variables
extern size_t nthreads;
extern volatile bool running;
extern int verbose;
extern uint64_t txn_flags;
extern double scale_factor;
extern uint64_t runtime;
extern uint64_t ops_per_worker;
extern int run_mode;
extern int enable_parallel_loading;
extern int pin_cpus;
extern int slow_exit;
extern int retry_aborted_transaction;
extern int no_reset_counters;
extern int backoff_aborted_transaction;

class scoped_db_thread_ctx {
public:
  scoped_db_thread_ctx(const scoped_db_thread_ctx &) = delete;
  scoped_db_thread_ctx(scoped_db_thread_ctx &&) = delete;
  scoped_db_thread_ctx &operator=(const scoped_db_thread_ctx &) = delete;

  scoped_db_thread_ctx(ndb_wrapper *db, bool loader)
    : db(db)
  {
    db->thread_init(loader);
  }
  ~scoped_db_thread_ctx()
  {
    db->thread_end();
  }
private:
  ndb_wrapper *const db;
};

class bench_loader : public ndb_thread {
public:
  bench_loader(unsigned long seed, ndb_wrapper *db,
               const std::map<std::string, ndb_ordered_index *> &open_tables)
    : r(seed), db(db), open_tables(open_tables), b(0),
      loader_id_(-1), loader_tids_(nullptr)
  {
    txn_obj_buf.reserve(str_arena::MinStrReserveLength);
    txn_obj_buf.resize(db->sizeof_txn_object(txn_flags));
  }

  inline void
  set_barrier(spin_barrier &b)
  {
    ALWAYS_ASSERT(!this->b);
    this->b = &b;
  }

  inline void pin_to_numa_nodes() {
#ifdef USE_BCC
    struct bitmask *bm = numa_allocate_nodemask();
    for (size_t isock = 0; isock < ::allocator::nsocks_used; ++isock) {
      numa_bitmask_setbit(bm,
          numa_node_of_cpu(::allocator::coreids_mapping[isock][0]));
    }
    ALWAYS_ASSERT(!numa_run_on_node_mask(bm));
    ALWAYS_ASSERT(!sched_yield());
    numa_free_nodemask(bm);
#else
    struct bitmask *bm = numa_allocate_nodemask();
    size_t nsocks = nthreads / 8;
    if (nthreads % 8 > 0)
      ++nsocks;
    for (size_t isock = 0; isock < nsocks; ++isock) {
      numa_bitmask_setbit(bm,
          numa_node_of_cpu(::allocator::coreids_mapping[isock][0]));
    }
    ALWAYS_ASSERT(!numa_run_on_node_mask(bm));
    ALWAYS_ASSERT(!sched_yield());
    numa_free_nodemask(bm);
#endif
  }

  virtual void
  run()
  {
    { 
      scoped_rcu_region r; // register this thread in rcu region
    }
    pin_to_numa_nodes();
    ALWAYS_ASSERT(b);
    b->count_down();
    b->wait_for();
    scoped_db_thread_ctx ctx(db, true);
    load();
  }

  void set_bcc_info(int loader_id, uint64_t *loader_tids)
  {
    loader_id_ = loader_id;
    loader_tids_ = loader_tids;
  }

protected:
  inline void *txn_buf() { return (void *) txn_obj_buf.data(); }

  virtual void load() = 0;

  util::fast_random r;
  ndb_wrapper *const db;
  std::map<std::string, ndb_ordered_index *> open_tables;
  spin_barrier *b;
  std::string txn_obj_buf;
  str_arena arena;

  int loader_id_;
  uint64_t *loader_tids_;
};

class bench_worker : public ndb_thread {
public:
  bench_worker(unsigned int worker_id,
               bool set_core_id,
               unsigned long seed, ndb_wrapper *db,
               const std::map<std::string, ndb_ordered_index *> &open_tables,
               spin_barrier *barrier_a, spin_barrier *barrier_b)
    : worker_id(worker_id), set_core_id(set_core_id),
      r(seed), db(db), open_tables(open_tables),
      barrier_a(barrier_a), barrier_b(barrier_b),
      // the ntxn_* numbers are per worker
      ntxn_commits(0), ntxn_aborts(0),
      latency_numer_us(0),
      backoff_shifts(0), // spin between [0, 2^backoff_shifts) times before retry
      size_delta(0),
      launched_threads_(0),
      nworkers_(0),
      //iworker_(0),
      //sh_conflict_(nullptr),
      tid_arrays_temp_(nullptr),
      last_tid_(0),
      sh_txn_ar_temp_(nullptr),
      worker_mem_(nullptr)
  {
    txn_obj_buf.reserve(str_arena::MinStrReserveLength);
    txn_obj_buf.resize(db->sizeof_txn_object(txn_flags));
  }

  virtual ~bench_worker() {}

  // returns [did_commit?, size_increase_bytes]
  typedef std::pair<bool, ssize_t> txn_result;
  typedef txn_result (*txn_fn_t)(bench_worker *);

  struct workload_desc {
    workload_desc() {}
    workload_desc(const std::string &name, double frequency, txn_fn_t fn)
      : name(name), frequency(frequency), fn(fn)
    {
      ALWAYS_ASSERT(frequency > 0.0);
      ALWAYS_ASSERT(frequency <= 1.0);
    }
    std::string name;
    double frequency;
    txn_fn_t fn;
  };
  typedef std::vector<workload_desc> workload_desc_vec;
  virtual workload_desc_vec get_workload() const = 0;

  virtual void run();

  inline size_t get_ntxn_commits() const { return ntxn_commits; }
  inline size_t get_ntxn_aborts() const { return ntxn_aborts; }

  inline uint64_t get_latency_numer_us() const { return latency_numer_us; }

  inline double
  get_avg_latency_us() const
  {
    return double(latency_numer_us) / double(ntxn_commits);
  }

  std::map<std::string, size_t> get_txn_counts() const;

  typedef ndb_wrapper::counter_map counter_map;
  typedef ndb_wrapper::txn_counter_map txn_counter_map;

#ifdef ENABLE_BENCH_TXN_COUNTERS
  inline txn_counter_map
  get_local_txn_counters() const
  {
    return local_txn_counters;
  }
#endif

  inline ssize_t get_size_delta() const { return size_delta; }

  void set_bcc_info(int launched_threads, int nworkers,
      uint64_t **temp_tid_arrays, uint64_t last_tid, shtxn **sh_txn_ar_temp)
  {
    launched_threads_ = launched_threads;
    nworkers_ = nworkers;
    //iworker_ = iworker;
    //sh_conflict_ = sh_conflict;
    tid_arrays_temp_ = temp_tid_arrays;
    last_tid_ = last_tid;
    sh_txn_ar_temp_ = sh_txn_ar_temp;
  }

  void bcc_txn_setup(void *txn) {
    int32_t windex = get_worker_index();
    uint64_t *start_tid = get_start_tid();

#ifdef BCC_READ_STARTTID_OPT
    if (worker_mem_->prev_tids_on_txn_validation_valid) {
      memcpy(start_tid, get_end_tid(), sizeof(uint64_t) * nworkers_);
      worker_mem_->prev_tids_on_txn_validation_valid = false;
    } else {
      read_worker_last_tids(start_tid);
    }
#else
    read_worker_last_tids(start_tid);
#endif
    uint64_t min_tid = ULONG_MAX;
    for (int i = 0; i < nworkers_; ++i) {
      if (i == windex)
        continue;
      if (start_tid[i] < min_tid)
        min_tid = start_tid[i];
    }

    shtxn *sh_txn = get_shtxn();
    sh_txn->release(ndb_wrapper::to_release_tid(min_tid));
    size_t rssize = db->size_of_hashed_read_set(txn);
    rssize = (rssize + 7UL) & (~7UL);   // 8-byte aligned
    auto allocret = get_shmm()->allocate(rssize);
    void *hv = db->initialize_hashed_read_set(txn, allocret.first);
    sh_txn->insert(
#if defined(BCC_READ_STARTTID_OPT) && defined(USE_BCC)
        worker_mem_->tid_arrays[::allocator::wid_to_sockid(windex)][::allocator::wid_to_idx_in_sock(windex)],
#else
        start_tid[windex],
#endif
        hv, allocret.first + rssize, allocret.second);

    db->set_txn_bcc_info(txn, nworkers_, windex, launched_threads_,
        start_tid, get_end_tid(),
        &worker_mem_->prev_tids_on_txn_validation_valid,
        get_shtxn_ar(), worker_mem_->tid_arrays, worker_mem_->unique_non_inserts);
  }

#ifdef BCC_TXN_STATS
  void bcc_txn_stats_setup(void *txn, const char *txn_name) {
    std::map<std::string, bcc_txn_stats>::iterator it = stats.find(txn_name);
    if (unlikely(it == stats.end())) {
      it = stats.insert(it, {std::string(txn_name), bcc_txn_stats(txn_name)});
    }
    db->set_stats_obj(txn, &it->second);
  }
#endif

  uint64_t *get_start_tid()
  {
    ALWAYS_ASSERT(worker_mem_);
    return worker_mem_->tids_on_txn_startup;
  }

  uint64_t *get_end_tid()
  {
    ALWAYS_ASSERT(worker_mem_);
    return worker_mem_->tids_on_txn_validation;
  }

  shtxn *get_shtxn()
  {
    ALWAYS_ASSERT(worker_mem_);
    return &worker_mem_->txn;
  }

  shtxn **get_shtxn_ar()
  {
    ALWAYS_ASSERT(worker_mem_);
    return worker_mem_->txns_ar;
  }

  shmm *get_shmm()
  {
    ALWAYS_ASSERT(worker_mem_);
    return &worker_mem_->mm;
  }

  int get_launched_threads()
  {
    return launched_threads_;
  }

  int get_nworkers()
  {
    return nworkers_;
  }

  int get_worker_index()
  {
    return worker_id % coreid::static_cpus_online;
//    return iworker_;
  }

//  uint64_t *get_tid_address()
//  {
//    return sh_tid_;
//  }

//  int64_t *get_sh_conflict()
//  {
//    return sh_conflict_;
//  }

  virtual void join() OVERRIDE
  {
    ndb_thread::join();
#ifdef BCC_STATS
    std::cerr << "Worker " << worker_id % coreid::static_cpus_online
              << " STATS:" << std::endl;
#endif
#ifdef USE_BCC
    get_shmm()->print_stats();
    get_shtxn()->print_stats();
#endif
#ifdef BCC_TXN_STATS
    for (auto &s : stats) {
      s.second.print();
    }
#endif
  }

protected:
  void read_worker_last_tids(uint64_t *tids_out) {
#ifdef USE_BCC
    uint64_t tids_read[BCC_MAX_NTHREADS];
    size_t ntids_read = 0;
    for (size_t isock = 0; isock < ::allocator::nsocks_used; ++isock) {
      memcpy(&tids_read[ntids_read], worker_mem_->tid_arrays[isock],
          sizeof(uint64_t) * ::allocator::nworkers_on_sock[isock]);
      ntids_read += ::allocator::nworkers_on_sock[isock];
    }
    ALWAYS_ASSERT(static_cast<int>(ntids_read) == nworkers_);
    // Shuffle tids_read into the order of worker index
    for (int wid = 0; wid < nworkers_; ++wid)
      tids_out[wid] = tids_read[::allocator::wid_to_idx_in_tids_read[wid]];
#endif
  }

#ifdef USE_BCC
  void InitializeWorkerMemoryFirstHalf()
  {
    const size_t wid = worker_id % coreid::static_cpus_online;
    worker_mem_ = reinterpret_cast<worker_numa_memory_layout *>(
        ::allocator::GetWorkerMemStart(wid));
    memset(worker_mem_->tids_on_txn_startup, 0,
        sizeof(worker_mem_->tids_on_txn_startup));
    memset(worker_mem_->tids_on_txn_validation, 0,
        sizeof(worker_mem_->tids_on_txn_validation));
    worker_mem_->prev_tids_on_txn_validation_valid = false;
    if (::allocator::is_wid_smallest_on_sock(wid)) {
      for (int iw = 0; iw < 8; ++iw)
        worker_mem_->local_tid_array[iw] = last_tid_;
      tid_arrays_temp_[::allocator::wid_to_sockid(wid)] = worker_mem_->local_tid_array;
    }
    worker_mem_->mm.init(wid, worker_mem_->mem,
        ::allocator::GetMemPerWorker());
    worker_mem_->txn.init(&worker_mem_->mm);
    sh_txn_ar_temp_[wid] = &worker_mem_->txn;
  }

  void InitializeWorkerMemorySecondHalf()
  {
    memcpy(worker_mem_->tid_arrays, tid_arrays_temp_, sizeof(uint64_t *) * 8);
    memcpy(worker_mem_->txns_ar, sh_txn_ar_temp_, sizeof(shtxn *) * nworkers_);
  }
#endif

  void bcc_pin_cpu(size_t core)
  {
    auto node = numa_node_of_cpu(core);
    // pin to node
    ALWAYS_ASSERT(!numa_run_on_node(node));
    // is numa_run_on_node() guaranteed to take effect immediately?
    ALWAYS_ASSERT(!sched_yield());
  }

#ifndef USE_BCC
  size_t wid_to_core(size_t wid) {
    size_t nsocks = nworkers_ / 8;
    if (nworkers_ % 8 > 0)
      ++nsocks;
    return ::allocator::coreids_mapping[wid % nsocks][wid / nsocks];
  }
#endif

// Make sure that BCC and Silo workers are both numa-pinned in the same way.
#define ALWAYS_PIN_CPU_ON_RUN_SETUP

  virtual void on_run_setup()
  {
#if defined(ALWAYS_PIN_CPU_ON_RUN_SETUP) || defined(USE_BCC)
    const size_t wid = worker_id % coreid::static_cpus_online; //wid_to_coreid(iworker_);
    if (!pin_cpus) {
#ifdef USE_BCC
      bcc_pin_cpu(::allocator::wid_to_coreid(wid));
#else
      bcc_pin_cpu(wid_to_core(wid));
#endif
    }
#endif
#ifdef USE_BCC
    ::allocator::FaultWorkerMemory(wid, ::allocator::wid_to_coreid(wid)/*, iworker_*/);
    InitializeWorkerMemoryFirstHalf();
#endif
  }

  void on_run_start()
  {
#ifdef USE_BCC
    InitializeWorkerMemorySecondHalf();
#endif
  }

  inline void *txn_buf() { return (void *) txn_obj_buf.data(); }

  unsigned int worker_id;
  bool set_core_id;
  util::fast_random r;
  ndb_wrapper *const db;
  std::map<std::string, ndb_ordered_index *> open_tables;
  spin_barrier *const barrier_a;
  spin_barrier *const barrier_b;

private:
  size_t ntxn_commits;
  size_t ntxn_aborts;
  uint64_t latency_numer_us;
  unsigned backoff_shifts;

protected:

#ifdef ENABLE_BENCH_TXN_COUNTERS
  txn_counter_map local_txn_counters;
  void measure_txn_counters(void *txn, const char *txn_name);
#else
  inline ALWAYS_INLINE void measure_txn_counters(void *txn, const char *txn_name) {}
#endif

  std::vector<size_t> txn_counts; // breakdown of txns
  ssize_t size_delta; // how many logical bytes (of values) did the worker add to the DB

  std::string txn_obj_buf;
  str_arena arena;

  int launched_threads_;
  int nworkers_;  // Total number of workers launched
  //int iworker_;   // The index of this worker
  //int64_t *sh_conflict_;
  uint64_t **tid_arrays_temp_;
  uint64_t last_tid_;  // Used to initialize local_tid_array
  shtxn **sh_txn_ar_temp_;
  worker_numa_memory_layout *worker_mem_;
#ifdef BCC_TXN_STATS
  std::map<std::string, bcc_txn_stats> stats;
#endif
};

class bench_runner {
public:
  bench_runner(const bench_runner &) = delete;
  bench_runner(bench_runner &&) = delete;
  bench_runner &operator=(const bench_runner &) = delete;

  bench_runner(ndb_wrapper *db, size_t nthreads)
    : db(db), barrier_a(nthreads), barrier_b(1), barrier_c(nthreads) {}
  virtual ~bench_runner() {}
  void run();
protected:
  virtual std::vector<bench_loader*> make_loaders() = 0;

  virtual std::vector<bench_worker*> make_workers() = 0;

  static void print_bcc_configs() {
    std::cerr << "BCC Configs:" << std::endl;
    std::cerr << "SHMM_MEMPOOL_SIZE_0: " << SHMM_MEMPOOL_SIZE_0 << std::endl;
    std::cerr << "SHMM_MEMPOOL_SIZE_1: " << SHMM_MEMPOOL_SIZE_1 << std::endl;
    std::cerr << "SHMM_MEMPOOL_SIZE_2: " << SHMM_MEMPOOL_SIZE_2 << std::endl;
    std::cerr << "MAX_TXN_ELEMS: " << MAX_TXN_ELEMS << std::endl;
    std::cerr << "READTID_OPT: "
#ifdef BCC_READ_STARTTID_OPT
              << true
#else
              << false
#endif
              << std::endl;
    std::cerr << "BCC_WR_LEFT_CONFLICT_CHECK_FORWARD: "
#ifdef BCC_WR_LEFT_CONFLICT_CHECK_FORWARD
              << true
#else
              << false
#endif
              << std::endl;
  }
  ndb_wrapper *const db;
  std::map<std::string, ndb_ordered_index *> open_tables;

  spin_barrier barrier_a;
  spin_barrier barrier_b;
  spin_barrier barrier_c;
};

class limit_callback : public ndb_ordered_index::scan_callback {
public:
  limit_callback(ssize_t limit = -1)
    : limit(limit), n(0)
  {
    ALWAYS_ASSERT(limit == -1 || limit > 0);
  }

  virtual bool invoke(
      const char *keyp, size_t keylen,
      const std::string &value)
  {
    INVARIANT(limit == -1 || n < size_t(limit));
    values.emplace_back(std::string(keyp, keylen), value);
    return (limit == -1) || (++n < size_t(limit));
  }

  typedef std::pair<std::string, std::string> kv_pair;
  std::vector<kv_pair> values;

  const ssize_t limit;
private:
  size_t n;
};


class latest_key_callback : public ndb_ordered_index::scan_callback {
public:
  latest_key_callback(std::string &k, ssize_t limit = -1)
    : limit(limit), n(0), k(&k)
  {
    ALWAYS_ASSERT(limit == -1 || limit > 0);
  }

  virtual bool invoke(
      const char *keyp, size_t keylen,
      const std::string &value)
  {
    INVARIANT(limit == -1 || n < size_t(limit));
    k->assign(keyp, keylen);
    ++n;
    return (limit == -1) || (n < size_t(limit));
  }

  inline size_t size() const { return n; }
  inline std::string &kstr() { return *k; }

private:
  ssize_t limit;
  size_t n;
  std::string *k;
};

// explicitly copies keys, because btree::search_range_call() interally
// re-uses a single string to pass keys (so using standard string assignment
// will force a re-allocation b/c of shared ref-counting)
//
// this isn't done for values, because each value has a distinct string from
// the string allocator, so there are no mutations while holding > 1 ref-count
template <size_t N>
class static_limit_callback : public ndb_ordered_index::scan_callback {
public:
  // XXX: push ignore_key into lower layer
  static_limit_callback(str_arena *arena, bool ignore_key)
    : n(0), arena(arena), ignore_key(ignore_key)
  {
    static_assert(N > 0, "xx");
  }

  virtual bool invoke(
      const char *keyp, size_t keylen,
      const std::string &value)
  {
    INVARIANT(n < N);
    INVARIANT(arena->manages(&value));
    if (ignore_key) {
      values.emplace_back(nullptr, &value);
    } else {
      std::string * const s_px = arena->next();
      INVARIANT(s_px && s_px->empty());
      s_px->assign(keyp, keylen);
      values.emplace_back(s_px, &value);
    }
    return ++n < N;
  }

  inline size_t
  size() const
  {
    return values.size();
  }

  typedef std::pair<const std::string *, const std::string *> kv_pair;
  typename util::vec<kv_pair, N>::type values;

private:
  size_t n;
  str_arena *arena;
  bool ignore_key;
};

#endif /* _NDB_BENCH_H_ */

