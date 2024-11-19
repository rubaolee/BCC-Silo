#ifndef _NDB_WRAPPER_H_
#define _NDB_WRAPPER_H_

#include "../txn_btree.h"
#include "../txn_proto2_impl.h"
#include "../str_arena.h"
#include "../schema.h"
#include "../sh_macros.h"

//#include "../sh_ht.h"

  enum TxnProfileHint {
    HINT_DEFAULT,

    // ycsb profiles
    HINT_KV_GET_PUT, // KV workloads over a single key
    HINT_KV_RMW, // get/put over a single key
    HINT_KV_SCAN, // KV scan workloads (~100 keys)

    // tpcc profiles
    HINT_TPCC_NEW_ORDER,
    HINT_TPCC_PAYMENT,
    HINT_TPCC_DELIVERY,
    HINT_TPCC_ORDER_STATUS,
    HINT_TPCC_ORDER_STATUS_READ_ONLY,
    HINT_TPCC_STOCK_LEVEL,
    HINT_TPCC_STOCK_LEVEL_READ_ONLY,

    // tpcw profiles
    HINT_TPCW_DO_CART,
    HINT_TPCW_DO_BUY_CONFIRM,

    HINT_BCCTEST,
  };


namespace private_ {
  struct ndbtxn {
    TxnProfileHint hint;
    char buf[0];
  } PACKED;

  // XXX: doesn't check to make sure you are passing in an ndbtx
  // of the right hint
  template <typename Traits>
  struct cast_base {
    typedef transaction_proto2<Traits> type;
    inline ALWAYS_INLINE type *
    operator()(struct ndbtxn *p) const
    {
      return reinterpret_cast<type *>(&p->buf[0]);
    }
  };

  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_get_probe0, ndb_get_probe0_cg)
  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_put_probe0, ndb_put_probe0_cg)
  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_insert_probe0, ndb_insert_probe0_cg)
  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_scan_probe0, ndb_scan_probe0_cg)
  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_remove_probe0, ndb_remove_probe0_cg)
  STATIC_COUNTER_DECL(scopedperf::tsc_ctr, ndb_dtor_probe0, ndb_dtor_probe0_cg)
}

class ndb_ordered_index{
protected:
  typedef private_::ndbtxn ndbtxn;
  template <typename Traits>
    using cast = private_::cast_base<Traits>;

public:
  ndb_ordered_index(const std::string &name, size_t value_size_hint, bool mostly_append);

  class scan_callback {
  public:
    ~scan_callback() {}
    // XXX(stephentu): key is passed as (const char *, size_t) pair
    // because it really should be the string_type of the underlying
    // tree, but since abstract_ordered_index is not templated we can't
    // really do better than this for now
    //
    // we keep value as std::string b/c we have more control over how those
    // strings are generated
    virtual bool invoke(const char *keyp, size_t keylen,
                        const std::string &value){return true;}
  };


  bool get(
      void *txn,
      const std::string &key,
      std::string &value, size_t max_bytes_read/*, int table*/);
  const char * put(
      void *txn,
      const std::string &key,
      const std::string &value/*, int table*/);
  const char * put(
      void *txn,
      std::string &&key,
      std::string &&value/*, int table*/);
  const char *
  insert(void *txn,
         const std::string &key,
         const std::string &value/*, int table*/);
  const char *
  insert(void *txn,
         std::string &&key,
         std::string &&value/*, int table*/);

  void scan(
      void *txn,
      const std::string &start_key,
      const std::string *end_key,
      scan_callback &callback,
      str_arena *arena/*, int table*/);

  void fullscan(void *txn, 
                const std::string &start_key,
                const std::string *end_key,
                my_scan_callback &callback, 
                str_arena *arena/*, int table*/);

  void rscan(
      void *txn,
      const std::string &start_key,
      const std::string *end_key,
      scan_callback &callback,
      str_arena *arena/*, int table*/);
  void remove(
      void *txn,
      const std::string &key);
  void remove(
      void *txn,
      std::string &&key);

  size_t size() const;
  std::map<std::string, uint64_t> clear();
private:
  std::string name;
  txn_btree<transaction_proto2> btr;
};


class ndb_wrapper{
protected:
  typedef private_::ndbtxn ndbtxn;
  template <typename Traits>
    using cast = private_::cast_base<Traits>;

public:

  ndb_wrapper(
      const std::vector<std::string> &logfiles,
      const std::vector<std::vector<unsigned>> &assignments_given,
      bool call_fsync,
      bool use_compression,
      bool fake_writes, size_t &nthreads);

  ssize_t txn_max_batch_size() const{ return 1; }

  class ndb_abort_exception{};

  void
  do_txn_epoch_sync() const
  {
    txn_epoch_sync<transaction_proto2>::sync();
  }

  void
  do_txn_finish() const
  {
    txn_epoch_sync<transaction_proto2>::finish();
  }

  void
  thread_init(bool loader)
  {
    txn_epoch_sync<transaction_proto2>::thread_init(loader);
  }

  void
  thread_end()
  {
    txn_epoch_sync<transaction_proto2>::thread_end();
  }

  std::tuple<uint64_t, uint64_t, double>
  get_ntxn_persisted() const
  {
    return txn_epoch_sync<transaction_proto2>::compute_ntxn_persisted();
  }

  void
  reset_ntxn_persisted()
  {
    txn_epoch_sync<transaction_proto2>::reset_ntxn_persisted();
  }

  size_t
  sizeof_txn_object(uint64_t txn_flags) const;

  void *new_txn(
      uint64_t txn_flags,
      str_arena &arena,
      void *buf,
      TxnProfileHint hint = HINT_DEFAULT);

#ifdef BCC_TXN_STATS
  void set_stats_obj(void *txn, bcc_txn_stats *stats);
#endif
  void set_txn_bcc_info(void *txn, int32_t nthreads, int32_t index,
      int32_t launched, uint64_t *start, uint64_t*end, bool *prev_tids_valid,
      shtxn** sh_txn_ar, uint64_t **tid_arrays, const void **unique_non_inserts);
  size_t size_of_hashed_read_set(void *txn) const;
  void *initialize_hashed_read_set(void *txn, char *rsaddr);
  static inline uint64_t to_release_tid(uint64_t tid) {
    return tid & (~transaction_proto2_static::CoreMask);
  }

  bool commit_txn(void *txn, uint64_t * tid);

  void abort_txn(void *txn);
  void print_txn_debug(void *txn) const;
  std::map<std::string, uint64_t> get_txn_counters(void *txn) const;

  ndb_ordered_index * open_index(const std::string &name,
             size_t value_size_hint,
             bool mostly_append = false);

  void
  close_index(ndb_ordered_index *idx);

  typedef std::map<std::string, uint64_t> counter_map;
  typedef std::map<std::string, counter_map> txn_counter_map;
};
#endif /* _NDB_WRAPPER_H_ */
