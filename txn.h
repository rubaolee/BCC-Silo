#ifndef _NDB_TXN_H_
#define _NDB_TXN_H_

#include <malloc.h>
#include <stdint.h>
#include <sys/types.h>
#include <pthread.h>

#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include <limits>
#include <type_traits>
#include <tuple>

#include <unordered_map>

#include "amd64.h"
#include "sh_ht.h"
#include "sh_txn.h"
#include "btree_choice.h"
#include "core.h"
#include "counter.h"
#include "macros.h"
#include "varkey.h"
#include "util.h"
#include "rcu.h"
#include "thread.h"
#include "spinlock.h"
#include "small_unordered_map.h"
#include "static_unordered_map.h"
#include "static_vector.h"
#include "prefetch.h"
#include "tuple.h"
#include "scopedperf.hh"
#include "marked_ptr.h"
#include "ndb_type_traits.h"

enum conflict_enum {
  NOCONFLICT = 2001,
  WRCONFLICT = 2002,
  WWCONFLICT = 2003,
  RWCONFLICT = 2004,
};

// forward decl
template <template <typename> class Transaction, typename P>
  class base_txn_btree;

class transaction_unusable_exception {};
class transaction_read_only_exception {};

// XXX: hacky
extern std::string (*g_proto_version_str)(uint64_t v);

// base class with very simple definitions- nothing too exciting yet
class transaction_base {
  template <template <typename> class T, typename P>
    friend class base_txn_btree;
public:

  typedef dbtuple::tid_t tid_t;
  typedef dbtuple::size_type size_type;
  typedef dbtuple::string_type string_type;

  // TXN_EMBRYO - the transaction object has been allocated but has not
  // done any operations yet
  enum txn_state { TXN_EMBRYO, TXN_ACTIVE, TXN_COMMITED, TXN_ABRT, };

  enum {
    // use the low-level scan protocol for checking scan consistency,
    // instead of keeping track of absent ranges
    TXN_FLAG_LOW_LEVEL_SCAN = 0x1,

    // true to mark a read-only transaction- if a txn marked read-only
    // does a write, a transaction_read_only_exception is thrown and the
    // txn is aborted
    TXN_FLAG_READ_ONLY = 0x2,

    // XXX: more flags in the future, things like consistency levels
  };

#define ABORT_REASONS(x) \
    x(ABORT_REASON_NONE) \
    x(ABORT_REASON_USER) \
    x(ABORT_REASON_UNSTABLE_READ) \
    x(ABORT_REASON_FUTURE_TID_READ) \
    x(ABORT_REASON_NODE_SCAN_WRITE_VERSION_CHANGED) \
    x(ABORT_REASON_NODE_SCAN_READ_VERSION_CHANGED) \
    x(ABORT_REASON_WRITE_NODE_INTERFERENCE) \
    x(ABORT_REASON_INSERT_NODE_INTERFERENCE) \
    x(ABORT_REASON_READ_NODE_INTEREFERENCE) \
    x(ABORT_REASON_READ_ABSENCE_INTEREFERENCE)

  enum abort_reason {
#define ENUM_X(x) x,
    ABORT_REASONS(ENUM_X)
#undef ENUM_X
  };

  static const char *
  AbortReasonStr(abort_reason reason)
  {
    switch (reason) {
#define CASE_X(x) case x: return #x;
    ABORT_REASONS(CASE_X)
#undef CASE_X
    default:
      break;
    }
    ALWAYS_ASSERT(false);
    return 0;
  }

  transaction_base(uint64_t flags)
    : state(TXN_EMBRYO),
      reason(ABORT_REASON_NONE),
      flags(flags) {}

  transaction_base(const transaction_base &) = delete;
  transaction_base(transaction_base &&) = delete;
  transaction_base &operator=(const transaction_base &) = delete;

protected:
#define EVENT_COUNTER_DEF_X(x) \
  static event_counter g_ ## x ## _ctr;
  ABORT_REASONS(EVENT_COUNTER_DEF_X)
#undef EVENT_COUNTER_DEF_X

  static event_counter *
  AbortReasonCounter(abort_reason reason)
  {
    switch (reason) {
#define EVENT_COUNTER_CASE_X(x) case x: return &g_ ## x ## _ctr;
    ABORT_REASONS(EVENT_COUNTER_CASE_X)
#undef EVENT_COUNTER_CASE_X
    default:
      break;
    }
    ALWAYS_ASSERT(false);
    return 0;
  }

public:

  // only fires during invariant checking
  inline void
  ensure_active()
  {
    if (state == TXN_EMBRYO)
      state = TXN_ACTIVE;
    INVARIANT(state == TXN_ACTIVE);
  }

  inline uint64_t
  get_flags() const
  {
    return flags;
  }

protected:

  // the read set is a mapping from (tuple -> tid_read).
  // "write_set" is used to indicate if this read tuple
  // also belongs in the write set.
  struct read_record_t {
    constexpr read_record_t() : tuple(), t() {}
    constexpr read_record_t(const dbtuple *tuple, tid_t t)
      : tuple(tuple), t(t) {}
    inline const dbtuple *
    get_tuple() const
    {
      return tuple;
    }
    inline tid_t
    get_tid() const
    {
      return t;
    }
  private:
    const dbtuple *tuple;
    tid_t t;
  };

  friend std::ostream &
  operator<<(std::ostream &o, const read_record_t &r);

  // the write set is logically a mapping from (tuple -> value_to_write).
  struct write_record_t {
    enum {
      FLAGS_INSERT  = 0x1,
      FLAGS_DOWRITE = 0x1 << 1,
    };

    constexpr inline write_record_t()
      : tuple(), k(), r(), w(), btr()
    {}

    // all inputs are assumed to be stable
    inline write_record_t(dbtuple *tuple,
                          const string_type *k,
                          const void *r,
                          dbtuple::tuple_writer_t w,
                          concurrent_btree *btr,
                          bool insert)
      : tuple(tuple),
        k(k),
        r(r),
        w(w),
        btr(btr)
    {
      this->btr.set_flags(insert ? FLAGS_INSERT : 0);
    }
    inline dbtuple *
    get_tuple()
    {
      return tuple;
    }
    inline const dbtuple *
    get_tuple() const
    {
      return tuple;
    }
    inline bool
    is_insert() const
    {
      return btr.get_flags() & FLAGS_INSERT;
    }
    inline bool
    do_write() const
    {
      return btr.get_flags() & FLAGS_DOWRITE;
    }
    inline void
    set_do_write()
    {
      INVARIANT(!do_write());
      btr.or_flags(FLAGS_DOWRITE);
    }
    inline concurrent_btree *
    get_btree() const
    {
      return btr.get();
    }
    inline const string_type &
    get_key() const
    {
      return *k;
    }
    inline const void *
    get_value() const
    {
      return r;
    }
    inline dbtuple::tuple_writer_t
    get_writer() const
    {
      return w;
    }
  private:
    dbtuple *tuple;
    const string_type *k;
    const void *r;
    dbtuple::tuple_writer_t w;
    marked_ptr<concurrent_btree> btr; // first bit for inserted, 2nd for dowrite
  };

  friend std::ostream &
  operator<<(std::ostream &o, const write_record_t &r);

  // the absent set is a mapping from (btree_node -> version_number).
  struct absent_record_t { uint64_t version; };

  friend std::ostream &
  operator<<(std::ostream &o, const absent_record_t &r);

  struct dbtuple_write_info {
    enum {
      FLAGS_LOCKED = 0x1,
      FLAGS_INSERT = 0x1 << 1,
    };
    dbtuple_write_info() : tuple(), entry(nullptr), pos() {}
    dbtuple_write_info(dbtuple *tuple, write_record_t *entry,
                       bool is_insert, size_t pos)
      : tuple(tuple), entry(entry), pos(pos)
    {
      if (is_insert)
        this->tuple.set_flags(FLAGS_LOCKED | FLAGS_INSERT);
    }
    // XXX: for searching only
    explicit dbtuple_write_info(const dbtuple *tuple)
      : tuple(const_cast<dbtuple *>(tuple)), entry(), pos() {}
    inline dbtuple *
    get_tuple()
    {
      return tuple.get();
    }
    inline const dbtuple *
    get_tuple() const
    {
      return tuple.get();
    }
    inline ALWAYS_INLINE void
    mark_locked()
    {
      INVARIANT(!is_locked());
      tuple.or_flags(FLAGS_LOCKED);
      INVARIANT(is_locked());
    }
    inline ALWAYS_INLINE bool
    is_locked() const
    {
      return tuple.get_flags() & FLAGS_LOCKED;
    }
    inline ALWAYS_INLINE bool
    is_insert() const
    {
      return tuple.get_flags() & FLAGS_INSERT;
    }
    inline ALWAYS_INLINE
    bool operator<(const dbtuple_write_info &o) const
    {
      // the unique key is [tuple, !is_insert, pos]
      return tuple < o.tuple ||
             (tuple == o.tuple && !is_insert() < !o.is_insert()) ||
             (tuple == o.tuple && !is_insert() == !o.is_insert() && pos < o.pos);
    }
    marked_ptr<dbtuple> tuple;
    write_record_t *entry;
    size_t pos;
  };


  static event_counter g_evt_read_logical_deleted_node_search;
  static event_counter g_evt_read_logical_deleted_node_scan;
  static event_counter g_evt_dbtuple_write_search_failed;
  static event_counter g_evt_dbtuple_write_insert_failed;

  static event_counter evt_local_search_lookups;
  static event_counter evt_local_search_write_set_hits;
  static event_counter evt_dbtuple_latest_replacement;

  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe0, g_txn_commit_probe0_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe1, g_txn_commit_probe1_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe2, g_txn_commit_probe2_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe3, g_txn_commit_probe3_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe4, g_txn_commit_probe4_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe5, g_txn_commit_probe5_cg);
  CLASS_STATIC_COUNTER_DECL(scopedperf::tsc_ctr, g_txn_commit_probe6, g_txn_commit_probe6_cg);

  txn_state state;
  abort_reason reason;
  const uint64_t flags;
};

// type specializations
namespace private_ {
  template <>
  struct is_trivially_destructible<transaction_base::read_record_t> {
    static const bool value = true;
  };

  template <>
  struct is_trivially_destructible<transaction_base::write_record_t> {
    static const bool value = true;
  };

  template <>
  struct is_trivially_destructible<transaction_base::absent_record_t> {
    static const bool value = true;
  };

  template <>
  struct is_trivially_destructible<transaction_base::dbtuple_write_info> {
    static const bool value = true;
  };
}

inline ALWAYS_INLINE std::ostream &
operator<<(std::ostream &o, const transaction_base::read_record_t &r)
{
  //o << "[tuple=" << util::hexify(r.get_tuple())
  o << "[tuple=" << *r.get_tuple()
    << ", tid_read=" << g_proto_version_str(r.get_tid())
    << "]";
  return o;
}

inline ALWAYS_INLINE std::ostream &
operator<<(
    std::ostream &o,
    const transaction_base::write_record_t &r)
{
  o << "[tuple=" << r.get_tuple()
    << ", key=" << util::hexify(r.get_key())
    << ", value=" << util::hexify(r.get_value())
    << ", insert=" << r.is_insert()
    << ", do_write=" << r.do_write()
    << ", btree=" << r.get_btree()
    << "]";
  return o;
}

inline ALWAYS_INLINE std::ostream &
operator<<(std::ostream &o, const transaction_base::absent_record_t &r)
{
  o << "[v=" << r.version << "]";
  return o;
}

struct default_transaction_traits {
  static const size_t read_set_expected_size = SMALL_SIZE_MAP;
  static const size_t absent_set_expected_size = EXTRA_SMALL_SIZE_MAP;
  static const size_t write_set_expected_size = SMALL_SIZE_MAP;
  static const bool stable_input_memory = false;
  static const bool hard_expected_sizes = false; // true if the expected sizes are hard maximums
  static const bool read_own_writes = true; // if we read a key which we previous put(), are we guaranteed
                                            // to read our latest (uncommited) values? this comes at a
                                            // performance penality [you should not need this behavior to
                                            // write txns, since you *know* the values you inserted]
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = SMALL_SIZE_MAP;
  static const size_t hash_buckets_expected_size = SMALL_SIZE_MAP;
  static const bool use_hashed_tuple_set = false;
#endif

  typedef util::default_string_allocator StringAllocator;
};

struct default_stable_transaction_traits : public default_transaction_traits {
  static const bool stable_input_memory = true;
};

#ifdef BCC_TXN_STATS
struct bcc_txn_stats {
  bcc_txn_stats(const char *name) :
    txn_name(name), ntxns_silo_committed(0), us_txns_silo_committed(0.0),
    ntxns_bcc_committed(0), us_txns_bcc_committed(0.0),
    us_before_bcc_checks_committed(0.0), us_in_bcc_checks_committed(0.0),
    ntxns_conflict_aborted(0), us_txns_conflict_aborted(0.0),
    us_before_bcc_checks_aborted(0.0), us_in_bcc_checks_aborted(0.0),
    ntxns_wr_left_conflict(0), ntxns_ww_left_conflict(0),
    ntxns_rw_left_conflict(0),
//    us_probe(0.0),
    ntxns_probed(0),
    ntxns_other_aborted(0), us_txns_other_aborted(0.0)
//    , n_write_sets(0), total_ratios_unique_non_inserts(0.0)
  {}
  bcc_txn_stats() : bcc_txn_stats("TXN") {}

  void add_txn_silo_committed(double us) {
    ++ntxns_silo_committed;
    us_txns_silo_committed += us;
  }

  void add_txn_bcc_committed(double us) {
    ++ntxns_bcc_committed;
    us_txns_bcc_committed += us;
  }

  void add_before_bcc_checks_committed(double us) {
    us_before_bcc_checks_committed += us;
  }

  void add_in_bcc_checks_committed(double us) {
    us_in_bcc_checks_committed += us;
  }

  void add_txn_conflict_aborted(double us) {
    ++ntxns_conflict_aborted;
    us_txns_conflict_aborted += us;
  }

  void add_before_bcc_checks_aborted(double us) {
    us_before_bcc_checks_aborted += us;
  }

  void add_in_bcc_checks_aborted(double us) {
    us_in_bcc_checks_aborted += us;
  }

  void add_txn_wr_left_conflict() {
    ++ntxns_wr_left_conflict;
  }

  void add_txn_ww_left_conflict() {
    ++ntxns_ww_left_conflict;
  }

  void add_txn_rw_left_conflict() {
    ++ntxns_rw_left_conflict;
  }

  void add_probe(/*double us*/) {
//    us_probe += us;
    ++ntxns_probed;
  }

  void add_txn_other_aborted(double us) {
    ++ntxns_other_aborted;
    us_txns_other_aborted += us;
  }

//  void add_ratio_unique_non_insert(double ratio) {
//    ++n_write_sets;
//    total_ratios_unique_non_inserts += ratio;
//  }

  void print() {
    std::cerr << txn_name << " STATS: "
              << "us_per_silo_committed(" << us_per_silo_committed()
              << "/" << ntxns_silo_committed << ") "
#ifdef USE_BCC
              << "us_per_bcc_committed(" << us_per_bcc_committed()
              << "/" << ntxns_bcc_committed << ") "
              << "us_before_bcc_checks_per_bcc_committed(" << us_before_bcc_checks_per_bcc_committed()
              << "/" << ntxns_bcc_committed << ") "
              << "us_in_bcc_checks_per_bcc_committed(" << us_in_bcc_checks_per_bcc_committed()
              << "/" << ntxns_bcc_committed << ") "
#endif
              << "us_per_conflict_aborted(" << us_per_conflict_aborted()
              << "/" << ntxns_conflict_aborted << ") "
#ifdef USE_BCC
              << "us_before_bcc_checks_per_conflict_aborted(" << us_before_bcc_checks_per_conflict_aborted()
              << "/" << ntxns_conflict_aborted << ") "
              << "us_in_bcc_checks_per_conflict_aborted(" << us_in_bcc_checks_per_conflict_aborted()
              << "/" << ntxns_conflict_aborted << ") "
              << "ntxns_wr_left_conflict(" << ntxns_wr_left_conflict << ") "
              << "ntxns_ww_left_conflict(" << ntxns_ww_left_conflict << ") "
              << "ntxns_rw_left_conflict(" << ntxns_rw_left_conflict << ") "
#endif
              << "us_per_other_aborted(" << us_per_other_aborted()
              << "/" << ntxns_other_aborted << ") "
//              << "ratio_per_unique_non_insert(" << ratio_per_unique_non_insert()
//              << "/" << n_write_sets << ") "
              << "ntxns_probed(" << ntxns_probed << ") "
              << std::endl;
  }

  double us_per_silo_committed() const {
    return (ntxns_silo_committed > 0) ?
        (us_txns_silo_committed / ntxns_silo_committed) : 0.0;
  }

  double us_per_bcc_committed() const {
    return (ntxns_bcc_committed > 0) ?
        (us_txns_bcc_committed / ntxns_bcc_committed) : 0.0;
  }

  double us_before_bcc_checks_per_bcc_committed() const {
    return (ntxns_bcc_committed > 0) ?
        (us_before_bcc_checks_committed / ntxns_bcc_committed) : 0.0;
  }

  double us_in_bcc_checks_per_bcc_committed() const {
    return (ntxns_bcc_committed > 0) ?
        (us_in_bcc_checks_committed / ntxns_bcc_committed) : 0.0;
  }

  double us_per_conflict_aborted() const {
    return (ntxns_conflict_aborted > 0) ?
        (us_txns_conflict_aborted / ntxns_conflict_aborted) : 0.0;
  }

  double us_before_bcc_checks_per_conflict_aborted() const {
    return (ntxns_conflict_aborted > 0) ?
        (us_before_bcc_checks_aborted / ntxns_conflict_aborted) : 0.0;
  }

  double us_in_bcc_checks_per_conflict_aborted() const {
    return (ntxns_conflict_aborted > 0) ?
        (us_in_bcc_checks_aborted / ntxns_conflict_aborted) : 0.0;
  }

//  double us_probe_per_txn() const {
//    return (ntxns_probed > 0) ?
//        (us_probe / ntxns_probed) : 0.0;
//  }

  double us_per_other_aborted() const {
    return (ntxns_other_aborted > 0) ?
        (us_txns_other_aborted / ntxns_other_aborted) : 0.0;
  }

//  double ratio_per_unique_non_insert() const {
//    return (n_write_sets > 0) ?
//        (total_ratios_unique_non_inserts / n_write_sets) : 0.0;
//  }

  const char *txn_name;
  // For txns that should have been committed even with vanilla Silo.
  uint64_t ntxns_silo_committed;
  double us_txns_silo_committed;
  // For txns that are aborted with vanilla Silo but committed with BCC.
  uint64_t ntxns_bcc_committed;
  double us_txns_bcc_committed;
  // Time before the beginning of BCC checks for txns aborted with vanilla
  // Silo but committed with BCC. This info is useful for analyzing the
  // tradeoffs between aborting a txn and doing BCC checks to save a
  // silo-aborted txn.
  // NOTE: The number of such txns is recorded in ntxns_bcc_committed.
  double us_before_bcc_checks_committed;
  // Time spent in BCC checks for txns aborted with vanilla Silo but
  // committed with BCC.
  // NOTE: The number of such txns is recorded in ntxns_bcc_committed.
  double us_in_bcc_checks_committed;
  // For txns aborted due to right- and left-edge conflict checks
  // (NO need to distinguish Silo/BCC in this case).
  uint64_t ntxns_conflict_aborted;
  double us_txns_conflict_aborted;
  // Time before the beginning of BCC checks for txns aborted even after
  // right- and left-edge conflict checks.
  // NOTE: The number of such txns is recorded in ntxns_conflict_aborted.
  double us_before_bcc_checks_aborted;
  // Time spent in BCC checks for txns aborted even after right- and
  // left-edge conflict checks.
  // NOTE: The number of such txns is recorded in ntxns_conflict_aborted.
  double us_in_bcc_checks_aborted;
  // Number of txns conflict_aborted due to wr left_check failures.
  uint64_t ntxns_wr_left_conflict;
  // Number of txns conflict_aborted due to ww left_check failures.
  uint64_t ntxns_ww_left_conflict;
  // Number of txns conflict_aborted due to rw left_check failures.
  uint64_t ntxns_rw_left_conflict;
  // NOTE: ntxns_wr_left_conflict + ntxns_ww_left_conflict +
  // ntxns_rw_left_conflict == ntxns_conflict_aborted
//  // Temporary: for probing.
//  double us_probe;
  uint64_t ntxns_probed;
  // For txns aborted due to all other reasons, such as WRITE_INTERFERENCE.
  uint64_t ntxns_other_aborted;
  double us_txns_other_aborted;
  // For computing percentage of unique non-insert tuples in write sets
//  uint64_t n_write_sets;
//  double total_ratios_unique_non_inserts;
};
#endif

template <template <typename> class Protocol, typename Traits>
class transaction : public transaction_base {
  // XXX: weaker than necessary
  template <template <typename> class, typename>
    friend class base_txn_btree;
  friend Protocol<Traits>;

public:

  // KeyWriter is expected to implement:
  // [1-arg constructor]
  //   KeyWriter(const Key *)
  // [fully materialize]
  //   template <typename StringAllocator>
  //   const std::string * fully_materialize(bool, StringAllocator &)

  // ValueWriter is expected to implement:
  // [1-arg constructor]
  //   ValueWriter(const Value *, ValueInfo)
  // [compute new size from old value]
  //   size_t compute_needed(const uint8_t *, size_t)
  // [fully materialize]
  //   template <typename StringAllocator>
  //   const std::string * fully_materialize(bool, StringAllocator &)
  // [perform write]
  //   void operator()(uint8_t *, size_t)
  //
  // ValueWriter does not have to be move/copy constructable. The value passed
  // into the ValueWriter constructor is guaranteed to be valid throughout the
  // lifetime of a ValueWriter instance.

  // KeyReader Interface
  //
  // KeyReader is a simple transformation from (const std::string &) => const Key &.
  // The input is guaranteed to be stable, so it has a simple interface:
  //
  //   const Key &operator()(const std::string &)
  //
  // The KeyReader is expect to preserve the following property: After a call
  // to operator(), but before the next, the returned value is guaranteed to be
  // valid and remain stable.

  // ValueReader Interface
  //
  // ValueReader is a more complex transformation from (const uint8_t *, size_t) => Value &.
  // The input is not guaranteed to be stable, so it has a more complex interface:
  //
  //   template <typename StringAllocator>
  //   bool operator()(const uint8_t *, size_t, StringAllocator &)
  //
  // This interface returns false if there was not enough buffer space to
  // finish the read, true otherwise.  Note that this interface returning true
  // does NOT mean that a read was stable, but it just means there were enough
  // bytes in the buffer to perform the tentative read.
  //
  // Note that ValueReader also exposes a dup interface
  //
  //   template <typename StringAllocator>
  //   void dup(const Value &, StringAllocator &)
  //
  // ValueReader also exposes a means to fetch results:
  //
  //   Value &results()
  //
  // The ValueReader is expected to preserve the following property: After a
  // call to operator(), if it returns true, then the value returned from
  // results() should remain valid and stable until the next call to
  // operator().

  //typedef typename P::Key key_type;
  //typedef typename P::Value value_type;
  //typedef typename P::ValueInfo value_info_type;

  //typedef typename P::KeyWriter key_writer_type;
  //typedef typename P::ValueWriter value_writer_type;

  //typedef typename P::KeyReader key_reader_type;
  //typedef typename P::SingleValueReader single_value_reader_type;
  //typedef typename P::ValueReader value_reader_type;

  typedef Traits traits_type;
  typedef typename Traits::StringAllocator string_allocator_type;

protected:
  // data structures

  inline ALWAYS_INLINE Protocol<Traits> *
  cast()
  {
    return static_cast<Protocol<Traits> *>(this);
  }

  inline ALWAYS_INLINE const Protocol<Traits> *
  cast() const
  {
    return static_cast<const Protocol<Traits> *>(this);
  }

  // XXX: we have baked in b-tree into the protocol- other indexes are possible
  // but we would need to abstract it away. we don't bother for now.

#ifdef USE_SMALL_CONTAINER_OPT
  // XXX: use parameterized typedef to avoid duplication

  // small types
  typedef small_vector<
    read_record_t,
    traits_type::read_set_expected_size> read_set_map_small;
  typedef small_vector<
    write_record_t,
    traits_type::write_set_expected_size> write_set_map_small;
  typedef small_unordered_map<
    const typename concurrent_btree::node_opaque_t *, absent_record_t,
    traits_type::absent_set_expected_size> absent_set_map_small;

  // static types
  typedef static_vector<
    read_record_t,
    traits_type::read_set_expected_size> read_set_map_static;
  typedef static_vector<
    write_record_t,
    traits_type::write_set_expected_size> write_set_map_static;
  typedef static_unordered_map<
    const typename concurrent_btree::node_opaque_t *, absent_record_t,
    traits_type::absent_set_expected_size> absent_set_map_static;

  // helper types for log writing
  typedef small_vector<
    uint32_t,
    traits_type::write_set_expected_size> write_set_u32_vec_small;
  typedef static_vector<
    uint32_t,
    traits_type::write_set_expected_size> write_set_u32_vec_static;

  // use static types if the expected sizes are guarantees
  typedef
    typename std::conditional<
      traits_type::hard_expected_sizes,
      read_set_map_static, read_set_map_small>::type read_set_map;
#ifdef USE_BCC
  typedef
    typename std::conditional<
      traits_type::use_hashed_tuple_set,
      HashedTupleSet<traits_type::hashed_tuple_set_expected_size,
        traits_type::hash_buckets_expected_size>,
      VectorTupleSet<traits_type::hashed_tuple_set_expected_size>>::type hashed_tuple_set;
#endif
  typedef
    typename std::conditional<
      traits_type::hard_expected_sizes,
      write_set_map_static, write_set_map_small>::type write_set_map;
  typedef
    typename std::conditional<
      traits_type::hard_expected_sizes,
      absent_set_map_static, absent_set_map_small>::type absent_set_map;
  typedef
    typename std::conditional<
      traits_type::hard_expected_sizes,
      write_set_u32_vec_static, write_set_u32_vec_small>::type write_set_u32_vec;

#else
  typedef std::vector<read_record_t> read_set_map;
  typedef std::vector<write_record_t> write_set_map;
  typedef std::vector<absent_record_t> absent_set_map;
  typedef std::vector<uint32_t> write_set_u32_vec;
#endif

  template <typename T>
    using write_set_sized_vec =
      typename std::conditional<
        traits_type::hard_expected_sizes,
        static_vector<T, traits_type::write_set_expected_size>,
        typename util::vec<T, traits_type::write_set_expected_size>::type
      >::type;

  // small type
  typedef
    typename util::vec<
      dbtuple_write_info, traits_type::write_set_expected_size>::type
    dbtuple_write_info_vec_small;

  // static type
  typedef
    static_vector<
      dbtuple_write_info, traits_type::write_set_expected_size>
    dbtuple_write_info_vec_static;

  // chosen type
  typedef
    typename std::conditional<
      traits_type::hard_expected_sizes,
      dbtuple_write_info_vec_static, dbtuple_write_info_vec_small>::type
    dbtuple_write_info_vec;

  static inline bool
  sorted_dbtuples_contains(
      const dbtuple_write_info_vec &dbtuples,
      const dbtuple *tuple)
  {
    // XXX: skip binary search for small-sized dbtuples?
    return std::binary_search(
        dbtuples.begin(), dbtuples.end(),
        dbtuple_write_info(tuple),
        [](const dbtuple_write_info &lhs, const dbtuple_write_info &rhs)
          { return lhs.get_tuple() < rhs.get_tuple(); });
  }

public:

  inline transaction(uint64_t flags, string_allocator_type &sa);
  inline ~transaction();

  // returns TRUE on successful commit, FALSE on abort
  // if doThrow, signals success by returning true, and
  // failure by throwing an abort exception
  bool commit(bool doThrow = false, uint64_t * tid = nullptr);

  // abort() always succeeds
  inline void
  abort()
  {
    abort_impl(ABORT_REASON_USER);
  }

  void dump_debug_info() const;

#ifdef DIE_ON_ABORT
  void
  abort_trap(abort_reason reason)
  {
    AbortReasonCounter(reason)->inc();
    this->reason = reason; // for dump_debug_info() to see
    dump_debug_info();
    ::abort();
  }
#else
  inline ALWAYS_INLINE void
  abort_trap(abort_reason reason)
  {
    AbortReasonCounter(reason)->inc();
  }
#endif

  std::map<std::string, uint64_t> get_txn_counters() const;

  inline ALWAYS_INLINE bool
  is_snapshot() const
  {
    return get_flags() & TXN_FLAG_READ_ONLY;
  }

  // for debugging purposes only
  inline const read_set_map &
  get_read_set() const
  {
    return read_set;
  }

  inline const write_set_map &
  get_write_set() const
  {
    return write_set;
  }

  inline const absent_set_map &
  get_absent_set() const
  {
    return absent_set;
  }

protected:
  inline void abort_impl(abort_reason r);

  // assumes lock on marker is held on marker by caller, and marker is the
  // latest: removes marker from tree, and clears latest
  void cleanup_inserted_tuple_marker(
      dbtuple *marker, const std::string &key,
      concurrent_btree *btr);

  // low-level API for txn_btree

  // try to insert a new "tentative" tuple into the underlying
  // btree associated with the given context.
  //
  // if return.first is not null, then this function will
  //   1) mutate the transaction such that the absent_set is aware of any
  //      mutating changes made to the underlying btree.
  //   2) add the new tuple to the write_set
  //
  // if return.second is true, then this txn should abort, because a conflict
  // was detected w/ the absent_set.
  //
  // if return.first is not null, the returned tuple is locked()!
  //
  // if the return.first is null, then this function has no side effects.
  //
  // NOTE: !ret.first => !ret.second
  // NOTE: assumes key/value are stable
  std::pair< dbtuple *, bool >
  try_insert_new_tuple(
      concurrent_btree &btr,
      const std::string *key,
      const void *value,
      dbtuple::tuple_writer_t writer/*, int table = 0*/);

  // reads the contents of tuple into v
  // within this transaction context
  template <typename ValueReader>
  bool
  do_tuple_read(const dbtuple *tuple, ValueReader &value_reader/*, int table = 0*/);

  void
  do_node_read(const typename concurrent_btree::node_opaque_t *n, uint64_t version);

public:
  // expected public overrides

  /**
   * Can we overwrite prev with cur?
   */
  bool can_overwrite_record_tid(tid_t prev, tid_t cur) const;

  inline string_allocator_type &
  string_allocator()
  {
    return *sa;
  }

public:
  /*
   * Update the committed tid info for the thread.
   * Called when a transaction is going to be committed.
   * Implemented in transaction_proto2.
   * Epoch should make progress when txn aborts.
   */
#ifdef USE_BCC
  inline void set_sh_tid(uint64_t tid) {
    //return;
    if (!tid_arrays)
      return;
#ifdef BCC_READ_STARTTID_OPT
    uint64_t old_tid = tid_arrays[::allocator::wid_to_sockid(tindex)][::allocator::wid_to_idx_in_sock(tindex)];
#else
    uint64_t old_tid = start_tid[tindex];
#endif
    uint64_t new_tid;
    if (tid > old_tid)
      new_tid = tid;
    else
      new_tid = old_tid + (static_cast<uint64_t>(1) << (cast()->get_num_id_shift()));
    tid_arrays[::allocator::wid_to_sockid(tindex)][::allocator::wid_to_idx_in_sock(tindex)] = new_tid;
  }
#endif

#ifdef BCC_TXN_STATS
  void set_stats_obj(bcc_txn_stats *set_stats) {
    stats = set_stats;
  }
#endif

  inline void set_bcc_info(
      int32_t nworkers, int32_t windex, int32_t launched,
      uint64_t *start, uint64_t *end, bool *prev_tids_valid,
      shtxn **stxn_ar, uint64_t **the_tid_arrays,
      const void **unique_non_inserts) {
    nworkers_ = nworkers;
    tindex = windex;
    launched_threads = launched;

    start_tid = start;
    end_tid = end;
    tids_valid = prev_tids_valid;
    left_conflict = NOCONFLICT;
    right_conflict = NOCONFLICT;

    sh_txn_ar = stxn_ar;
    sh_txn = sh_txn_ar[tindex];
    tid_arrays = the_tid_arrays;

    buf_unique_non_inserts = unique_non_inserts;
}

  // NOTE: Should only be called ifdef USE_BCC.
  size_t size_of_hashed_read_set() const {
#ifdef USE_BCC
    return sizeof(*hashed_read_set);
#endif
    // This should not be reached.
    return 0;
  }

  void *initialize_hashed_read_set(char *rsaddr) {
    void *hv = nullptr;
#ifdef USE_BCC
    hashed_read_set = new(rsaddr) hashed_tuple_set();
    hv = static_cast<TupleSearchInterface *>(hashed_read_set);
#endif
    return hv;
  }

protected:
#ifdef BCC_WR_LEFT_CONFLICT_CHECK_FORWARD
  bool check_rw_right_conflict() {
    if (right_conflict != NOCONFLICT)
      return true;
    if (!read_set.empty()) {
      typename read_set_map::iterator it     = read_set.begin();
      typename read_set_map::iterator it_end = read_set.end();
      for (; it != it_end; ++it) {
        if (find_write_set2(it->get_tuple())) {
          if(it->get_tuple()->is_latest_version(it->get_tid()))
            continue;
        } else {
          if (it->get_tuple()->stable_is_latest_version(it->get_tid()))
            continue;
        }
        right_conflict = RWCONFLICT;
        return true;
      }
    }
    return false;
  }

  bool check_wr_left_conflict2(tid_t rtid) {
    if (rtid == dbtuple::MAX_TID)
      return false;
    uint64_t cid = cast()->get_core_id(rtid);  // core id
    uint64_t widx = cid % coreid::static_cpus_online;  // worker index
    if (cid <= static_cast<uint64_t>(launched_threads) ||
        widx == static_cast<uint64_t>(tindex))
      return false;
#ifdef BCC_CHECK
    ALWAYS_ASSERT(widx >= 0 && widx < static_cast<uint64_t>(nworkers_));
#endif
    if (rtid > start_tid[widx])
      return true;
    return false;
  }
#endif

  bool check_wr_left_conflict() {
    ALWAYS_ASSERT(start_tid);
    typename read_set_map::iterator it     = read_set.begin();
    typename read_set_map::iterator it_end = read_set.end();
    for (; it != it_end; ++it) {
      tid_t rtid = it->get_tid();
      if (rtid == dbtuple::MAX_TID)
        continue;
      uint64_t cid = cast()->get_core_id(rtid);  // core id
      uint64_t widx = cid % coreid::static_cpus_online;  // worker index
      if (cid <= static_cast<uint64_t>(launched_threads) ||
          widx == static_cast<uint64_t>(tindex))
        continue;
#ifdef BCC_CHECK
      ALWAYS_ASSERT(widx >= 0 && widx < static_cast<uint64_t>(nworkers_));
#endif
      if (rtid > start_tid[widx])
        return true;
    }
    return false;
  }

  // XXX: use sorted dbwrite_set to avoid checking for repeated writes
  bool check_ww_left_conflict() {
    ALWAYS_ASSERT(start_tid);
    typename write_set_map::iterator it     = write_set.begin();
    typename write_set_map::iterator it_end = write_set.end();
    for (; it != it_end; ++it) {
      if (it->is_insert())
        continue;
      // The tuple should have been locked, so we can directly read its tid.
#ifdef BCC_CHECK
      ALWAYS_ASSERT(it->get_tuple());
#endif
      tid_t wtid = it->get_tuple()->version;
      if (wtid == dbtuple::MAX_TID)
        continue;
      uint64_t cid = cast()->get_core_id(wtid);  // core id
      uint64_t widx = cid % coreid::static_cpus_online;  // worker index
      if (cid <= static_cast<uint64_t>(launched_threads) ||
          widx == static_cast<uint64_t>(tindex))
        continue;
#ifdef BCC_CHECK
      ALWAYS_ASSERT(widx >= 0 && widx < static_cast<uint64_t>(nworkers_));
#endif
      if (wtid > start_tid[widx])
        return true;
    }
    return false;
  }

  // This incurs the highest overhead in BCC checking because
  // it does lots of inter-worker memory reading.
  // XXX: use sorted dbwrite_set to avoid checking for repeated writes
  bool check_rw_left_conflict() {
    ALWAYS_ASSERT(tid_arrays && end_tid && start_tid && sh_txn_ar);
    read_worker_last_tids(end_tid);
    is_end_tid_populated = true;
    *tids_valid = true;
    // For SSE comparison
    if (n_unique_non_inserts & 1)
      buf_unique_non_inserts[n_unique_non_inserts] = nullptr;
    for (int32_t i = 0; i < nworkers_; ++i) {
      if (i == tindex)
        continue;
//      typename write_set_map::iterator it = write_set.begin();
//      typename write_set_map::iterator it_end = write_set.end();
//      for (; it != it_end; ++it) {
//        if (it->is_insert())
//          continue;
//#ifdef BCC_CHECK
//        ALWAYS_ASSERT(it->get_tuple());
//#endif
      if (sh_txn_ar[i]->check_rw_conflict(start_tid[i],
          end_tid[i], buf_unique_non_inserts, n_unique_non_inserts/*it->get_tuple()*/))
        return true;
//      }
    }
    return false;
  }

  void read_worker_last_tids(uint64_t *tids_out) {
#ifdef USE_BCC
    uint64_t tids_read[BCC_MAX_NTHREADS];
    size_t ntids_read = 0;
    for (size_t isock = 0; isock < ::allocator::nsocks_used; ++isock) {
      memcpy(&tids_read[ntids_read], tid_arrays[isock],
          sizeof(uint64_t) * ::allocator::nworkers_on_sock[isock]);
      ntids_read += ::allocator::nworkers_on_sock[isock];
    }
    ALWAYS_ASSERT(static_cast<int32_t>(ntids_read) == nworkers_);
    // Shuffle tids_read into the order of worker index
    for (int32_t wid = 0; wid < nworkers_; ++wid)
      tids_out[wid] = tids_read[::allocator::wid_to_idx_in_tids_read[wid]];
#endif
  }

  uint64_t get_release_tid() {
    ALWAYS_ASSERT(tid_arrays && end_tid);
    uint64_t release_tid = 0;
    if (!is_end_tid_populated) {
      read_worker_last_tids(end_tid);
      is_end_tid_populated = true;
      *tids_valid = true;
    }
    for (int i = 0; i < nworkers_; ++i) {
      if (i == tindex)
        continue;
      if (end_tid[i] > release_tid)
        release_tid = end_tid[i];
    }
    return release_tid & (~(cast()->get_core_mask()));
  }

protected:
  // expected protected overrides

  /**
   * create a new, unique TID for a txn. at the point which gen_commit_tid(),
   * it still has not been decided whether or not this txn will commit
   * successfully
   */
  tid_t gen_commit_tid(const dbtuple_write_info_vec &write_tuples);
  tid_t gen_abort_tid();

  bool can_read_tid(tid_t t) const;

  // For GC handlers- note that on_dbtuple_spill() is called
  // with the lock on ln held, to simplify GC code
  //
  // Is also called within an RCU read region
  void on_dbtuple_spill(dbtuple *tuple_ahead, dbtuple *tuple);

  // Called when the latest value written to ln is an empty
  // (delete) marker. The protocol can then decide how to schedule
  // the logical node for actual deletion
  void on_logical_delete(dbtuple *tuple, const std::string &key, concurrent_btree *btr);

  // if gen_commit_tid() is called, then on_tid_finish() will be called
  // with the commit tid. before on_tid_finish() is called, state is updated
  // with the resolution (commited, aborted) of this txn
  void on_tid_finish(tid_t commit_tid);

  void on_post_rcu_region_completion();

protected:
  inline void clear();

  // SLOW accessor methods- used for invariant checking

  typename read_set_map::iterator
  find_read_set(const dbtuple *tuple)
  {
    // linear scan- returns the *first* entry found
    // (a tuple can exist in the read_set more than once)
    typename read_set_map::iterator it     = read_set.begin();
    typename read_set_map::iterator it_end = read_set.end();
    for (; it != it_end; ++it)
      if (it->get_tuple() == tuple)
        break;
    return it;
  }

  inline typename read_set_map::const_iterator
  find_read_set(const dbtuple *tuple) const
  {
    return const_cast<transaction *>(this)->find_read_set(tuple);
  }

  typename write_set_map::iterator
  find_write_set(dbtuple *tuple)
  {
    // linear scan- returns the *first* entry found
    // (a tuple can exist in the write_set more than once)
    typename write_set_map::iterator it     = write_set.begin();
    typename write_set_map::iterator it_end = write_set.end();
    for (; it != it_end; ++it)
      if (it->get_tuple() == tuple)
        break;
    return it;
  }

#ifdef BCC_WR_LEFT_CONFLICT_CHECK_FORWARD
  bool find_write_set2(const dbtuple *tuple)
  {
    // linear scan- returns the *first* entry found
    // (a tuple can exist in the write_set more than once)
    typename write_set_map::iterator it     = write_set.begin();
    typename write_set_map::iterator it_end = write_set.end();
    for (; it != it_end; ++it)
      if (it->get_tuple() == tuple)
        return true;
    return false;
  }
#endif

  inline typename write_set_map::const_iterator
  find_write_set(const dbtuple *tuple) const
  {
    return const_cast<transaction *>(this)->find_write_set(tuple);
  }

  inline bool
  handle_last_tuple_in_group(
      dbtuple_write_info &info, bool did_group_insert);

  read_set_map read_set;
  write_set_map write_set;
  absent_set_map absent_set;

  string_allocator_type *sa;

  shtxn *sh_txn;      // Local worker's txn manager
  shtxn **sh_txn_ar;  // Global txn manager array
  uint64_t **tid_arrays;   // Global per-sock tid arrays
  //int64_t * sh_conflict;
  uint64_t *start_tid;
  uint64_t *end_tid;
  bool *tids_valid;
  const void **buf_unique_non_inserts;
  int n_unique_non_inserts;
  bool is_end_tid_populated;
#ifdef USE_BCC
  hashed_tuple_set *hashed_read_set;
#endif
  int32_t nworkers_;          /* total number of worker threads in the system */
  int32_t tindex;             /* index of the current worker thread */
  int32_t launched_threads;
  int left_conflict;   // Indicates whether there is wr/ww/rw conflict in L->T
  int right_conflict;  // Indicates whether there is r/w conflict in T->R

#ifdef BCC_TXN_STATS
  double t_txn_start;
  double t_bcc_check_start;
  double t_bcc_check_end;
  bool probed;
  bcc_txn_stats *stats;
#endif

  unmanaged<scoped_rcu_region> rcu_guard_;
};

class transaction_abort_exception : public std::exception {
public:
  transaction_abort_exception(transaction_base::abort_reason r)
    : r(r) {}
  inline transaction_base::abort_reason
  get_reason() const
  {
    return r;
  }
  virtual const char *
  what() const throw()
  {
    return transaction_base::AbortReasonStr(r);
  }
private:
  transaction_base::abort_reason r;
};

// XXX(stephentu): stupid hacks
// XXX(stephentu): txn_epoch_sync is a misnomer
template <template <typename> class Transaction>
struct txn_epoch_sync {
  // block until the next epoch
  static inline void sync() {}
  // finish any async jobs
  static inline void finish() {}
  // run this code when a benchmark worker finishes
  static inline void thread_end() {}
  // how many txns have we persisted in total, from
  // the last reset invocation?
  static inline std::pair<uint64_t, double>
    compute_ntxn_persisted() { return {0, 0.0}; }
  // reset the persisted counters
  static inline void reset_ntxn_persisted() {}
};

#endif /* _NDB_TXN_H_ */
