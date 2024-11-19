#ifndef __BCC_TXN_H__
#define __BCC_TXN_H__

#include <climits>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <vector>

#include "macros.h"
#include "sh_ht.h"
#include "sh_macros.h"
#include "sh_mm.h"

#define NEXT_N_TXN_IDX(i, n)  (((i) + (n)) & (MAX_TXN_ELEMS - 1))
#define NEXT_TXN_IDX(i)       NEXT_N_TXN_IDX(i, 1)
#define PREV_TXN_IDX(i)       (((i) + MAX_TXN_ELEMS - 1) & (MAX_TXN_ELEMS - 1))

#ifdef BCC_SHTXN_STATS
struct shtxn_stats {
  void init() {
    max_list_len = 0;
    curr_list_len = 0;
    total_list_len = 0;
    ninserts = 0;
    max_scanned = 0;
    total_scanned = 0;
    nscans = 0;
#ifdef BCC_SHHT_STATS
    total_fullness_committed = 0.0;
    ntxns_committed = 0;
    total_fullness_aborted = 0.0;
    ntxns_aborted = 0;
#endif
  }

  void print() {
    std::cerr << "SHTXN STATS: "
              << "avg_list_len("
              << (ninserts > 0 ? (total_list_len / ninserts) : 0) << ") "
              << "max_list_len(" << max_list_len << ") "
              << "avg_rw_conflict_scanned("
              << (nscans > 0 ? (total_scanned / nscans) : 0) << "/"
              << nscans << ") "
              << "max_rw_conflict_scanned(" << max_scanned << "/"
              << nscans << ") "
#ifdef BCC_SHHT_STATS
              << "fullness_per_txn(" << fullness_per_txn()
              << "/" << (ntxns_committed + ntxns_aborted) << ") "
              << "fullness_per_committed(" << fullness_per_committed()
              << "/" << ntxns_committed << ") "
              << "fullness_per_aborted(" << fullness_per_aborted()
              << "/" << ntxns_aborted << ") "
#endif
              << std::endl;
  }

  void add_insert() {
    if (++curr_list_len > max_list_len)
      max_list_len = curr_list_len;
    total_list_len += curr_list_len;
    ++ninserts;
  }

  void add_release(size_t release_count) {
    curr_list_len -= release_count;
  }

  void add_scan(size_t scanned) {
    if (scanned > max_scanned)
      max_scanned = scanned;
    total_scanned += scanned;
    ++nscans;
  }

#ifdef BCC_SHHT_STATS
  void add_fullness_committed(double fullness) {
    total_fullness_committed += fullness;
    ++ntxns_committed;
  }

  void add_fullness_aborted(double fullness) {
    total_fullness_aborted += fullness;
    ++ntxns_aborted;
  }

  double fullness_per_txn() {
    double total_fullness = total_fullness_committed + total_fullness_aborted;
    size_t ntxns = ntxns_committed + ntxns_aborted;
    return ntxns > 0 ? (total_fullness / ntxns) : 0.0;
  }

  double fullness_per_committed() {
    return ntxns_committed > 0 ? (total_fullness_committed / ntxns_committed)
        : 0.0;
  }

  double fullness_per_aborted() {
    return ntxns_aborted > 0 ? (total_fullness_aborted / ntxns_aborted)
        : 0.0;
  }
#endif

  // Maximum length of txn list during the worker's lifetime.
  size_t max_list_len;
  // Current length of txn list; used only for computing max_list_len.
  // *** Do NOT print this field ***
  size_t curr_list_len;
  // For computing average length of txn list.
  size_t total_list_len;
  uint64_t ninserts;
  // For computing average/max scan length in check_rw_conflict returning false.
  size_t max_scanned;
  size_t total_scanned;
  uint64_t nscans;
#ifdef BCC_SHHT_STATS
  // Total fullness of all committed txns
  double total_fullness_committed;
  size_t ntxns_committed;
  // Total fullness of all aborted txns
  double total_fullness_aborted;
  size_t ntxns_aborted;
#endif
};
#endif

/*
 * Each worker maintains a list of processed transactions
 * to help other peer workers reduce false aborts.
 * NOTE: Do NOT create an object of this class on the stack; only on heap!
 */
class shtxn {
protected:
  struct txn_list_elem {
   // Last committed tid of the worker when this txn starts.
    uint64_t last_tid;
    // When all peer workers' last_tids become larger than release_tid,
    // this txn elem can be released.
    // This field will be set to 0 if this txn aborts in the end.
    uint64_t release_tid;
    // NOTE: hv may not equal the starting address of this txn's
    // underlying read set memory area maintained by sh_mm.
    const TupleSearchInterface *hv;
    // Offset of the right margin of this txn's shmm area relative to hv.
    // NOTE: This assumes that the size of each txn's hashed_tuple_set can
    // be expressed with an unsigned short integer.
    unsigned short end_addr_offset;
    // Which mempool
    char ipool;
  };

  // txn_elems[start_pos, end_pos) are occupied.
  struct txn_list_head {
    uint32_t start_pos;
    uint32_t end_pos;
  };

public:
  inline void init(shmm *shm) {
    memset(&txn_head, 0, sizeof(txn_head));
    sh_mm = shm;
#ifdef BCC_SHTXN_STATS
    stats.init();
#endif
  }

  // Add a new in-flight txn.
  inline void insert(uint64_t last_tid, const void *hv,
      const void *end_addr, int ipool) {
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    return;
#endif
#ifdef BCC_CHECK
    ALWAYS_ASSERT((is_empty() || is_release_tid_set()) && !is_full());
    ALWAYS_ASSERT(end_addr > hv);
#endif
    uint32_t pos = txn_head.end_pos;
    txn_elems[pos].last_tid = last_tid;
    txn_elems[pos].release_tid = ULONG_MAX;
    txn_elems[pos].hv = reinterpret_cast<const TupleSearchInterface *>(hv);
    txn_elems[pos].end_addr_offset = static_cast<unsigned short>(
        reinterpret_cast<uint64_t>(end_addr) - reinterpret_cast<uint64_t>(hv));
    txn_elems[pos].ipool = static_cast<char>(ipool);
    COMPILER_MEMORY_FENCE;
    txn_head.end_pos = NEXT_TXN_IDX(pos);
#ifdef BCC_SHTXN_STATS
    stats.add_insert();
#endif
  }

  /*
   * Set the release tid of the in-flight txn.
   */
  inline void set_release(uint64_t tid, bool aborted = false) {
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    return;
#endif
#ifdef BCC_CHECK
    ALWAYS_ASSERT(tid != ULONG_MAX && !is_release_tid_set());
#endif
    uint32_t pos = PREV_TXN_IDX(txn_head.end_pos);
    if (aborted) {
      //txn_elems[pos].last_tid = 0;
#ifdef BCC_SHHT_STATS
      stats.add_fullness_aborted(txn_elems[pos].hv->Fullness());
#endif
    }
#ifdef BCC_SHHT_STATS
    else {
      stats.add_fullness_committed(txn_elems[pos].hv->Fullness());
    }
#endif
    txn_elems[pos].release_tid = tid;
  }

  // Returns whether the release tid for the in-flight txn has been set.
  inline bool is_release_tid_set() const {
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    return false;
#endif
#ifdef BCC_CHECK
    ALWAYS_ASSERT(!is_empty());
#endif
    uint32_t pos = PREV_TXN_IDX(txn_head.end_pos);
    return txn_elems[pos].release_tid != ULONG_MAX;
  }

  inline bool check_rw_conflict(
      uint64_t start_tid, uint64_t end_tid, const void **tuples, int ntuples/*const void *tuple*/) {
#ifdef BCC_CHECK
    if (start_tid > end_tid)
      std::cerr << start_tid << " " << end_tid << std::endl;
    ALWAYS_ASSERT(start_tid <= end_tid);
#endif
#ifdef BCC_SHTXN_STATS
    size_t scanned = 0;
#endif
    // Take a view of current start_pos and end_pos.
    uint32_t prev_start_pos = PREV_TXN_IDX(txn_head.start_pos);
    uint32_t prev_end_pos = PREV_TXN_IDX(txn_head.end_pos);
    for (uint32_t k = prev_end_pos; k != prev_start_pos; k = PREV_TXN_IDX(k)) {
#ifdef BCC_SHTXN_STATS
      ++scanned;
#endif
      if (txn_elems[k].release_tid == 0 || txn_elems[k].last_tid > end_tid)
        continue;
      if (txn_elems[k].last_tid < start_tid)
        break;
      if (txn_elems[k].hv->SearchTuple(tuples, ntuples/*tuple*/))
        return true;
    }
#ifdef BCC_SHTXN_STATS
    stats.add_scan(scanned);
#endif
    return false;
  }

  /*
   * Release memory for txns whose read sets no longer need to be kept.
   * Txns whose release_tids are less than @tid can be released.
   */
  inline void release(uint64_t tid) {
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    return;
#endif
    size_t release_count = 0;
    const char *end_addr[SHMM_MEMPOOLS] = {nullptr};
    for (uint32_t i = txn_head.start_pos; i != txn_head.end_pos;
        i = NEXT_TXN_IDX(i)) {
      if (txn_elems[i].release_tid < tid) {
        release_count++;
#ifdef BCC_CHECK
        ALWAYS_ASSERT(txn_elems[i].ipool >= 0 &&
            txn_elems[i].ipool < SHMM_MEMPOOLS);
#endif
        end_addr[static_cast<int>(txn_elems[i].ipool)] =
            reinterpret_cast<const char *>(txn_elems[i].hv) +
            txn_elems[i].end_addr_offset;
      } else {
        break;
      }
    }
    if (release_count > 0) {
      txn_head.start_pos = NEXT_N_TXN_IDX(txn_head.start_pos, release_count);
      COMPILER_MEMORY_FENCE;
      for (int ipool = 0; ipool < SHMM_MEMPOOLS; ++ipool) {
        if (end_addr[ipool])
          sh_mm->release(ipool, end_addr[ipool]);
      }
#ifdef BCC_SHTXN_STATS
      stats.add_release(release_count);
#endif
    }
  }

  inline bool is_empty() const {
    return txn_head.end_pos == txn_head.start_pos;
  }

  inline bool is_full() const {
    return NEXT_TXN_IDX(txn_head.end_pos) == txn_head.start_pos;
  }

  void print_stats() {
#ifdef BCC_SHTXN_STATS
    stats.print();
#endif
  }

private:
  // A circular list of finished/in-flight txns.
  txn_list_head txn_head;
  txn_list_elem txn_elems[MAX_TXN_ELEMS];
  // The underlying memory manager.
  shmm *sh_mm;
#ifdef BCC_SHTXN_STATS
  shtxn_stats stats;
#endif
};

#endif
