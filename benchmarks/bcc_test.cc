#include <iostream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>

#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#include <sys/sysinfo.h>

#include "../macros.h"
#include "../varkey.h"
#include "../thread.h"
#include "../util.h"
#include "../spinbarrier.h"
#include "../base_txn_btree.h"
#include "../schema.h"

#include "../record/encoder.h"
#include "../record/inline_str.h"
#include "../sh_ht.h"
#include "../sh_txn.h"
#include "../sh_macros.h"

#include "bccbench.h"
#include "bccndb_wrapper_impl.h"

static uint64_t update_num = 1;
static uint64_t update_scan_num = 5;
static uint64_t hnum = 1;
static double hperc = 0;

using namespace std;
using namespace util;

static size_t nkeys;
class test_worker : public bench_worker {
public:
  test_worker(unsigned int worker_id,
               unsigned long seed, ndb_wrapper *db,
               const map<string, ndb_ordered_index *> &open_tables,
               spin_barrier *barrier_a, spin_barrier *barrier_b,
               uint64_t id)
    : bench_worker(worker_id, true, seed, db, open_tables, barrier_a, barrier_b),
      tbl(open_tables.at("testscan")), id(id), seed(seed) { }

  txn_result txn_update_query() {
    scoped_str_arena s_arena(arena);
    void *txn = db->new_txn(txn_flags, arena, txn_buf(), HINT_BCCTEST);
#ifdef USE_BCC
    bcc_txn_setup(txn);
#endif
#ifdef BCC_TXN_STATS
    bcc_txn_stats_setup(txn, "bcc_test");
#endif
    uint64_t start_pos = 500;
    try {
      int tindex = worker_id % coreid::static_cpus_online;

      if (tindex != 0) {
        for (uint32_t i = 0; i < update_scan_num; i++) {
            if(i ==0 && (r.next()%100)/100.0 < hperc){
                k1.t1 = 1;
            }else{
                k1.t1 = start_pos + i;
            }
            tbl->get(txn, Encode(k1), obj_key2, -1);
        }

        start_pos += update_scan_num;

        for (uint32_t i = 0; i < update_num; i++) {
          k1.t1 = 100 + tindex * update_num + i;
          v1.t2 = 1;
          tbl->put(txn, Encode(k1), Encode(v1));
        }

      } else {

        for (uint32_t i=0;i<hnum;i++){
          k1.t1 = i+1;
          v1.t2 = 1;
          tbl->put(txn,Encode(k1),Encode(v1));
        }
      }
      measure_txn_counters(txn,"bcc_test");

      if (likely(db->commit_txn(txn, NULL))) {
        txn_result ret(true, 10);
        return ret;
      } else {
        /* this means txn is aborted and no exception is thrown */
        db->abort_txn(txn);
        return txn_result(false, 0);
      }
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      db->abort_txn(txn);
    }
    return txn_result(false, 0);
  }

  static txn_result TxnUpdateQuery(bench_worker *w) {
    return static_cast<test_worker *>(w)->txn_update_query();
  }

  virtual workload_desc_vec get_workload() const {
    workload_desc_vec w;
    w.push_back(workload_desc("TestUpdateQuery", 1.0, TxnUpdateQuery));
    return w;
  }

protected:
  virtual void
  on_run_setup() OVERRIDE
  {
    const size_t a = worker_id % coreid::num_cpus_online();
    const size_t b = a % nthreads;
    if (!pin_cpus) {
      //return;
      goto finish;
    }
    rcu::s_instance.pin_current_thread(b);
    rcu::s_instance.fault_region();
  finish:
    bench_worker::on_run_setup();
  }


private:
    ndb_ordered_index *tbl;
    uint64_t id;
    unsigned long seed;
    vector<std::string> x;
    testscan::key k1,k2;
    testscan::value v1;
    string obj_key1, obj_key2;
};

class test_table_loader : public bench_loader {
public:
  test_table_loader(unsigned long seed,
                     ndb_wrapper *db,
                     const map<string, ndb_ordered_index *> &open_tables)
    : bench_loader(seed, db, open_tables) {}

protected:
  virtual void load() {
    ndb_ordered_index *tbl = open_tables.at("testscan");
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    try {
      // load
      const size_t batchsize = (db->txn_max_batch_size() == -1) ?
        10000 : db->txn_max_batch_size();
      ALWAYS_ASSERT(batchsize > 0);
      const size_t nbatches = nkeys / batchsize;
      for (size_t id = 0; id < nthreads; id++) {
        if (nbatches == 0) {
          void *txn = db->new_txn(txn_flags, arena, txn_buf());
          for (size_t j = 0; j < nkeys; j++) {
            testscan::key kk;
            testscan::value vv;
            kk.t1 = j;
            vv.t2 = j;
            tbl->insert(txn, Encode(kk),Encode(vv));
          }
          if (verbose)
            cerr << "batch 1/1 done" << endl;

#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn,&tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif

        } else {
          for (size_t i = 0; i < nbatches; i++) {
            size_t keyend = (i == nbatches - 1) ? nkeys : (i + 1) * batchsize;
            void *txn = db->new_txn(txn_flags, arena, txn_buf());
            for (size_t j = i * batchsize; j < keyend; j++) {
                testscan::key kk ;
                testscan::value vv ;
                kk.t1 = j;
                vv.t2 = j;
              tbl->insert(txn, Encode(kk), Encode(vv));
            }
            if (verbose)
              cerr << "batch " << (i + 1) << "/" << nbatches << " done" << endl;

#ifdef USE_BCC
            ALWAYS_ASSERT(db->commit_txn(txn,&tid));
            if (tid > latest_tid)
              latest_tid = tid;
#else
            ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif

          }
        }
      }
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }

#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose)
      cerr << "[INFO] finished loading table" << endl;
  }
};

class test_bench_runner : public bench_runner {
public:
  test_bench_runner(ndb_wrapper *db, size_t nthreads): bench_runner(db,nthreads){
    open_tables["testscan"] = db->open_index("testscan", 10);
  }

protected:
  virtual vector<bench_loader *>make_loaders() OVERRIDE {
    vector<bench_loader *> ret;
    ret.push_back(new test_table_loader(0, db, open_tables));
    return ret;
  }

  virtual vector<bench_worker *>make_workers() OVERRIDE {
    const unsigned alignment = coreid::num_cpus_online();
    const int blockstart =
      coreid::allocate_contiguous_aligned_block(nthreads, alignment);
    ALWAYS_ASSERT(blockstart >= 0);
    ALWAYS_ASSERT((blockstart % alignment) == 0);
    fast_random r(8544290);
    vector<bench_worker *> ret;
    for (size_t i = 0; i < nthreads; i++) {
        ret.push_back(
          new test_worker(
            blockstart + i, r.next(), db, open_tables,
            &barrier_a, &barrier_b, i));
      }
    return ret;
  }
};

void bcctest(ndb_wrapper *db, size_t nthreads, int argc, char**argv){

  nkeys = size_t(scale_factor * 1000.0);
  ALWAYS_ASSERT(nkeys > 0);
  optind = 1;
  while (1) {
    static struct option long_options[] =
    {
      {"update-scan-num"            , required_argument , 0                             , 'r'} ,
      {"hot-percent"            , required_argument , 0                             , 'h'} ,
      {"hot-update-num"                         , required_argument , 0                  , 'w'} ,
      {"update-num"                         , required_argument , 0                     , 'u'} ,
      {0, 0, 0, 0}
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "r:", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0)
        break;
      abort();
      break;

    case 'h':
        hperc = atoi(optarg)/100.0;
        break;

    case 'r':
      update_scan_num = strtoul(optarg, NULL, 10);
      break;

    case 'w':
        hnum = atoi(optarg);
        break;

    case 'u':
      update_num = strtoul(optarg, NULL, 10);
      break;

    case '?':
      /* getopt_long already printed an error message. */
      exit(1);

    default:
      abort();
    }
  }

  test_bench_runner r(db, nthreads);
  r.run();
}
