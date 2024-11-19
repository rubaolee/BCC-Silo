#ifndef _NDB_WRAPPER_IMPL_H_
#define _NDB_WRAPPER_IMPL_H_

#include <stdint.h>
#include "bccndb_wrapper.h"
#include "../counter.h"
#include "../rcu.h"
#include "../varkey.h"
#include "../macros.h"
#include "../util.h"
#include "../scopedperf.hh"
#include "../txn.h"
//#include "../txn_proto1_impl.h"
#include "../txn_proto2_impl.h"
#include "../tuple.h"

struct hint_default_traits : public default_transaction_traits {
  typedef str_arena StringAllocator;
};

struct hint_bcctest_traits {
  static const size_t read_set_expected_size = 16;
  static const size_t write_set_expected_size = 16;
  static const size_t absent_set_expected_size = read_set_expected_size / 7 + 1;
  static const bool stable_input_memory = false;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 16;
  static const size_t hash_buckets_expected_size = 4;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};


// ycsb profiles

struct hint_kv_get_put_traits {
  static const size_t read_set_expected_size = 16;
  static const size_t write_set_expected_size = 16;
  static const size_t absent_set_expected_size = read_set_expected_size / 7 + 1;
  static const bool stable_input_memory = false;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 1;   // Determined by -n
  static const size_t hash_buckets_expected_size = 1;       // Determined by -n
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};

struct hint_kv_rmw_traits : public hint_kv_get_put_traits {};

struct hint_kv_scan_traits {
  static const size_t read_set_expected_size = 100;
  static const size_t write_set_expected_size = 1;
  static const size_t absent_set_expected_size = read_set_expected_size / 7 + 1;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 1;   // Determined by -n
  static const size_t hash_buckets_expected_size = 1;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};

// tpcc profiles

struct hint_read_only_traits {
  static const size_t read_set_expected_size = 1;
  static const size_t write_set_expected_size = 1;
  static const size_t absent_set_expected_size = 1;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = true;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 1;
  static const size_t hash_buckets_expected_size = 1;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_new_order_traits {
  static const size_t read_set_expected_size = 35;
  static const size_t write_set_expected_size = 35;
  static const size_t absent_set_expected_size = 1;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = true;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 35;
  static const size_t hash_buckets_expected_size = 4;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_payment_traits {
  static const size_t read_set_expected_size = 85;
  static const size_t write_set_expected_size = 10;
  static const size_t absent_set_expected_size = 15;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  // Upper bound: 150; but can be set smaller
  static const size_t hashed_tuple_set_expected_size = 85;
  static const size_t hash_buckets_expected_size = 16;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_delivery_traits {
  static const size_t read_set_expected_size = 175;
  static const size_t write_set_expected_size = 175;
  static const size_t absent_set_expected_size = 35;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 175;
  static const size_t hash_buckets_expected_size = 32;
  static const bool use_hashed_tuple_set = true;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_order_status_traits {
  static const size_t read_set_expected_size = 95;
  static const size_t write_set_expected_size = 1;
  static const size_t absent_set_expected_size = 25;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  // Upper bound 514; but can be set smaller
  static const size_t hashed_tuple_set_expected_size = 95;
  static const size_t hash_buckets_expected_size = 16;
  static const bool use_hashed_tuple_set = true;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_order_status_read_only_traits : public hint_read_only_traits {};

struct hint_tpcc_stock_level_traits {
  static const size_t read_set_expected_size = 500;
  static const size_t write_set_expected_size = 1;
  static const size_t absent_set_expected_size = 25;
  static const bool stable_input_memory = true;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = false;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 550;
  static const size_t hash_buckets_expected_size = 64;
  static const bool use_hashed_tuple_set = true;
#endif
  typedef str_arena StringAllocator;
};

struct hint_tpcc_stock_level_read_only_traits : public hint_read_only_traits {};


struct hint_tpcw_do_cart_traits {
  static const size_t read_set_expected_size = 35;
  static const size_t write_set_expected_size = 35;
  static const size_t absent_set_expected_size = 25;
  static const bool stable_input_memory = false;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = true;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 10;
  static const size_t hash_buckets_expected_size = 4;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};


struct hint_tpcw_do_buy_confirm_traits {
  static const size_t read_set_expected_size = 35;
  static const size_t write_set_expected_size = 35;
  static const size_t absent_set_expected_size = 25;
  static const bool stable_input_memory = false;
  static const bool hard_expected_sizes = false;
  static const bool read_own_writes = true;
#ifdef USE_BCC
  static const size_t hashed_tuple_set_expected_size = 4;
  static const size_t hash_buckets_expected_size = 4;
  static const bool use_hashed_tuple_set = false;
#endif
  typedef str_arena StringAllocator;
};



#define TXN_PROFILE_HINT_OP(x) \
  x(HINT_DEFAULT, hint_default_traits) \
  x(HINT_KV_GET_PUT, hint_kv_get_put_traits) \
  x(HINT_BCCTEST, hint_bcctest_traits) \
  x(HINT_KV_RMW, hint_kv_rmw_traits) \
  x(HINT_KV_SCAN, hint_kv_scan_traits) \
  x(HINT_TPCC_NEW_ORDER, hint_tpcc_new_order_traits) \
  x(HINT_TPCC_PAYMENT, hint_tpcc_payment_traits) \
  x(HINT_TPCC_DELIVERY, hint_tpcc_delivery_traits) \
  x(HINT_TPCC_ORDER_STATUS, hint_tpcc_order_status_traits) \
  x(HINT_TPCC_ORDER_STATUS_READ_ONLY, hint_tpcc_order_status_read_only_traits) \
  x(HINT_TPCC_STOCK_LEVEL, hint_tpcc_stock_level_traits) \
  x(HINT_TPCC_STOCK_LEVEL_READ_ONLY, hint_tpcc_stock_level_read_only_traits)    \
  x(HINT_TPCW_DO_CART, hint_tpcw_do_cart_traits)    \
  x(HINT_TPCW_DO_BUY_CONFIRM, hint_tpcw_do_buy_confirm_traits)

ndb_wrapper::ndb_wrapper(
    const std::vector<std::string> &logfiles,
    const std::vector<std::vector<unsigned>> &assignments_given,
    bool call_fsync,
    bool use_compression,
    bool fake_writes, size_t &nthreads)
{
  if (logfiles.empty())
    return;
  std::vector<std::vector<unsigned>> assignments_used;
  txn_logger::Init(
      nthreads, logfiles, assignments_given, &assignments_used,
      call_fsync,
      use_compression,
      fake_writes);
}

size_t
ndb_wrapper::sizeof_txn_object(uint64_t txn_flags) const
{
#define MY_OP_X(a, b) sizeof(typename cast< b >::type),
  const size_t xs[] = {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  };
#undef MY_OP_X
  size_t xmax = 0;
  for (size_t i = 0; i < ARRAY_NELEMS(xs); i++)
    xmax = std::max(xmax, xs[i]);
  return xmax;
}

void *
ndb_wrapper::new_txn(
    uint64_t txn_flags,
    str_arena &arena,
    void *buf,
    TxnProfileHint hint)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(buf);
  p->hint = hint;
#define MY_OP_X(a, b) \
  case a: \
    new (&p->buf[0]) typename cast< b >::type(txn_flags, arena); \
    return p;
  switch (hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
  return 0;
}

template <typename T>
static inline ALWAYS_INLINE void
Destroy(T *t)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_dtor_probe0_cg);
  t->~T();
}

bool
ndb_wrapper::commit_txn(void *txn, uint64_t *tid = nullptr)
{
  // WARN(BCC): t->commit must not throw exception in case of transaction
  // aborts because the overhead of C++ exception handling is too high
  // (tens to hundreds of microseconds). That would cause workers spending
  // most of their time on handling exceptions thrown by aborted transactions
  // when contention is high (and transaction abortion rate is high).
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try{
#define MY_OP_X(a, b) \
      case a: \
        { \
          auto t = cast< b >()(p); \
          const bool ret = t->commit(false, tid); \
          if (ret) Destroy(t); \
          return ret; \
        }
      switch (p->hint) {
        TXN_PROFILE_HINT_OP(MY_OP_X)
      default:
        ALWAYS_ASSERT(false);
      }
    #undef MY_OP_X
      return false;
  }catch (transaction_abort_exception &ex) {
    throw ndb_abort_exception();
  }
}

#ifdef BCC_TXN_STATS
void ndb_wrapper::set_stats_obj(void *txn, bcc_txn_stats *stats) {
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      t->set_stats_obj(stats); \
      return; \
    }
  switch (p->hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
}
#endif

void ndb_wrapper::set_txn_bcc_info(void *txn, int32_t nthreads, int32_t index,
    int32_t launched, uint64_t *start, uint64_t *end, bool *prev_tids_valid,
    shtxn** sh_txn_ar, uint64_t **tid_arrays, const void **unique_non_inserts) {
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      t->set_bcc_info(nthreads, index, launched, start, end, prev_tids_valid, sh_txn_ar, tid_arrays, unique_non_inserts); \
      return; \
    }
  switch (p->hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
}

size_t ndb_wrapper::size_of_hashed_read_set(void *txn) const {
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
case a: \
  { \
    auto t = cast< b >()(p); \
    size_t rssize = t->size_of_hashed_read_set(); \
    return rssize; \
  }
switch (p->hint) {
  TXN_PROFILE_HINT_OP(MY_OP_X)
default:
  ALWAYS_ASSERT(false);
}
#undef MY_OP_X
  return 0;
}

void *ndb_wrapper::initialize_hashed_read_set(void *txn, char *rsaddr) {
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
case a: \
  { \
    auto t = cast< b >()(p); \
    void *hv = t->initialize_hashed_read_set(rsaddr); \
    return hv; \
  }
switch (p->hint) {
  TXN_PROFILE_HINT_OP(MY_OP_X)
default:
  ALWAYS_ASSERT(false);
}
#undef MY_OP_X
  return nullptr;
}

void
ndb_wrapper::abort_txn(void *txn)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      t->abort(); \
      Destroy(t); \
      return; \
    }
  switch (p->hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
}

void
ndb_wrapper::print_txn_debug(void *txn) const
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      t->dump_debug_info(); \
      return; \
    }
  switch (p->hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
}

std::map<std::string, uint64_t>
ndb_wrapper::get_txn_counters(void *txn) const
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      return t->get_txn_counters(); \
    }
  switch (p->hint) {
    TXN_PROFILE_HINT_OP(MY_OP_X)
  default:
    ALWAYS_ASSERT(false);
  }
#undef MY_OP_X
  return std::map<std::string, uint64_t>();
}

ndb_ordered_index *
ndb_wrapper::open_index(const std::string &name, size_t value_size_hint, bool mostly_append)
{
  return new ndb_ordered_index(name, value_size_hint, mostly_append);
}

void
ndb_wrapper::close_index(ndb_ordered_index *idx)
{
  //delete idx;
}

ndb_ordered_index::ndb_ordered_index(
    const std::string &name, size_t value_size_hint, bool mostly_append)
  : name(name), btr(value_size_hint, mostly_append, name)
{
  // for debugging
  //std::cerr << name << " : btree= "
  //          << btr.get_underlying_btree()
  //          << std::endl;
}

bool
ndb_ordered_index::get(
    void *txn,
    const std::string &key,
    std::string &value, size_t max_bytes_read/*, int table = 0*/)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_get_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      if (!btr.search(*t, key, value, max_bytes_read)) \
        return false; \
      return true; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
    INVARIANT(!value.empty());
    return true;
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception(); //return false;
  }
}

// XXX: find way to remove code duplication below using C++ templates!

const char *
ndb_ordered_index::put(
    void *txn,
    const std::string &key,
    const std::string &value/*, int table = 0*/)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_put_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.put(*t, key, value); \
      return 0; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
  return 0;
}

const char *
ndb_ordered_index::put(
    void *txn,
    std::string &&key,
    std::string &&value/*, int table = 0*/)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.put(*t, std::move(key), std::move(value)); \
      return 0; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
  return 0;
}

const char *
ndb_ordered_index::insert(
    void *txn,
    const std::string &key,
    const std::string &value/*, int table = 0*/)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_insert_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.insert(*t, key, value); \
      return 0; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
  return 0;
}

const char *
ndb_ordered_index::insert(
    void *txn,
    std::string &&key,
    std::string &&value/*, int table = 0*/)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.insert(*t, std::move(key), std::move(value)); \
      return 0; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
  return 0;
}

class ndb_wrapper_search_range_callback : public txn_btree<transaction_proto2>::search_range_callback {
public:
  ndb_wrapper_search_range_callback(ndb_ordered_index::scan_callback &upcall)
    : upcall(&upcall) {}

  virtual bool
  invoke(const typename txn_btree<transaction_proto2>::keystring_type &k,
         const typename txn_btree<transaction_proto2>::string_type &v)
  {
    return upcall->invoke(k.data(), k.length(), v);
  }

private:
  ndb_ordered_index::scan_callback *upcall;
};

/*
void
ndb_ordered_index::fullscan(
    void *txn,
    const std::string &start_key,
    const std::string *end_key,
    my_scan_callback &callback,
    str_arena *arena) {

  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_scan_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.fullscan(*t, start_key, end_key, callback); \
      return; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    ;
  }
}
*/

void
ndb_ordered_index::scan(
    void *txn,
    const std::string &start_key,
    const std::string *end_key,
    scan_callback &callback,
    str_arena *arena/*, int table = 0*/)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_scan_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  ndb_wrapper_search_range_callback c(callback);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.search_range_call(*t, start_key, end_key, c); \
      return; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
}

void
ndb_ordered_index::rscan(
    void *txn,
    const std::string &start_key,
    const std::string *end_key,
    scan_callback &callback,
    str_arena *arena/*, int table = 0*/)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  ndb_wrapper_search_range_callback c(callback);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.rsearch_range_call(*t, start_key, end_key, c); \
      return; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
}

void
ndb_ordered_index::remove(void *txn, const std::string &key)
{
  PERF_DECL(static std::string probe1_name(std::string(__PRETTY_FUNCTION__) + std::string(":total:")));
  ANON_REGION(probe1_name.c_str(), &private_::ndb_remove_probe0_cg);
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.remove(*t, key); \
      return; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
}

void
ndb_ordered_index::remove(void *txn, std::string &&key)
{
  ndbtxn * const p = reinterpret_cast<ndbtxn *>(txn);
  try {
#define MY_OP_X(a, b) \
  case a: \
    { \
      auto t = cast< b >()(p); \
      btr.remove(*t, std::move(key)); \
      return; \
    }
    switch (p->hint) {
      TXN_PROFILE_HINT_OP(MY_OP_X)
    default:
      ALWAYS_ASSERT(false);
    }
#undef MY_OP_X
  } catch (transaction_abort_exception &ex) {
    throw ndb_wrapper::ndb_abort_exception();
  }
}

size_t
ndb_ordered_index::size() const
{
  return btr.size_estimate();
}

std::map<std::string, uint64_t>
ndb_ordered_index::clear()
{
#ifdef TXN_BTREE_DUMP_PURGE_STATS
  std::cerr << "purging txn index: " << name << std::endl;
#endif
  return btr.unsafe_purge(true);
}

#endif /* _NDB_WRAPPER_IMPL_H_ */
