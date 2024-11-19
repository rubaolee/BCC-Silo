/**
 * An implementation of TPC-W subset transactions
 */

#include <sys/time.h>
#include <string>
#include <ctype.h>
#include <stdlib.h>
#include <malloc.h>

#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include <set>
#include <vector>

#include "../txn.h"
#include "../macros.h"
#include "../scopedperf.hh"
#include "../spinlock.h"

#include "bccbench.h"
#include "tpcw.h"
#include "bccndb_wrapper_impl.h"

using namespace std;
using namespace util;

#define TPCW_TABLE_LIST(x) \
  x(address) \
  x(author)  \
  x(cc_xacts)   \
  x(country)  \
  x(customer) \
  x(item)   \
  x(order_line) \
  x(orders) \
  x(shopping_cart)  \
  x(shopping_cart_line) \
  x(next_id)


static string credit_cards[] = {"VISA", "MASTERCARD", "DISCOVER",
                          "AMEX", "DINERS"};
static int num_card_types = 5;
string ship_types[] = {"AIR", "UPS", "FEDEX", "SHIP", "COURIER", "MAIL"};
static int num_ship_types = 6;
static string status_types[] = {"PROCESSING", "SHIPPED", "PENDING",
                         "DENIED"};
static int num_status_types = 4;

static unsigned g_txn_workload_mix[] = {100,0};
static double gtheta = 0.8; 
static int do_buy = 1;
static uint32_t rperc = 100;

struct _dummy {}; // exists so we can inherit from it, so we can use a macro in
                  // an init list...

class tpcw_worker_mixin : private _dummy {

#define DEFN_TBL_INIT_X(name) \
  , tbl_ ## name ## _vec(partitions.at(#name))

public:
  tpcw_worker_mixin(const map<string, vector<ndb_ordered_index *>> &partitions) :
    _dummy() // so hacky...
    TPCW_TABLE_LIST(DEFN_TBL_INIT_X)
  {
  }

#undef DEFN_TBL_INIT_X

protected:

#define DEFN_TBL_ACCESSOR_X(name) \
private:  \
  vector<ndb_ordered_index *> tbl_ ## name ## _vec; \
protected: \
  inline ALWAYS_INLINE ndb_ordered_index * \
  tbl_ ## name (unsigned int wid) \
  { \
    return tbl_ ## name ## _vec[wid]; \
  }

  TPCW_TABLE_LIST(DEFN_TBL_ACCESSOR_X)

#undef DEFN_TBL_ACCESSOR_X


public:

  static constexpr inline ALWAYS_INLINE int 
  NUM_EBS(){
    return 100;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_IDS(){
    return 64;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_SCL(){
    return 1000;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_ITEMS(){
    return 10000;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_CUSTOMERS(){
    return NUM_EBS() * 2880;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_ADDRESSES(){
    return NUM_CUSTOMERS() * 2;
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_AUTHORS(){
    return (int) (.25*NUM_ITEMS());
  }

  static constexpr inline ALWAYS_INLINE int 
  NUM_ORDERS(){
    return (int) (0.9 * NUM_CUSTOMERS());
  }

  static inline uint32_t
  GetCurrentTimeMillis()
  {
    //struct timeval tv;
    //ALWAYS_ASSERT(gettimeofday(&tv, 0) == 0);
    //return tv.tv_sec * 1000;

    // XXX(stephentu): implement a scalable GetCurrentTimeMillis()
    // for now, we just give each core an increasing number

    static __thread uint32_t tl_hack = 0;
    return tl_hack++;
  }

  // utils for generating random #s and strings

  static inline ALWAYS_INLINE int
  CheckBetweenInclusive(int v, int lower, int upper)
  {
    INVARIANT(v >= lower);
    INVARIANT(v <= upper);
    return v;
  }

  static inline ALWAYS_INLINE int
  RandomNumber(fast_random &r, int min, int max)
  {
    return CheckBetweenInclusive((int) (r.next_uniform() * (max - min + 1) + min), min, max);
  }

  static inline ALWAYS_INLINE int
  NonUniformRandom(fast_random &r, int A, int C, int min, int max)
  {
    return (((RandomNumber(r, 0, A) | RandomNumber(r, min, max)) + C) % (max - min + 1)) + min;
  }


  static string NameTokens[];

  // following oltpbench, we really generate strings of len - 1...
  static inline string
  RandomStr(fast_random &r, uint len)
  {
    // this is a property of the oltpbench implementation...
    if (!len)
      return "";

    uint i = 0;
    string buf(len - 1, 0);
    while (i < (len - 1)) {
      const char c = (char) r.next_char();
      // XXX(stephentu): oltpbench uses java's Character.isLetter(), which
      // is a less restrictive filter than isalnum()
      if (!isalnum(c))
        continue;
      buf[i++] = c;
    }
    return buf;
  }

  // RandomNStr() actually produces a string of length len
  static inline string
  RandomNStr(fast_random &r, uint len)
  {
    const char base = '0';
    string buf(len, 0);
    for (uint i = 0; i < len; i++)
      buf[i] = (char)(base + (r.next() % 10));
    return buf;
  }
};

string tpcw_worker_mixin::NameTokens[] =
  {
    string("BAR"),
    string("OUGHT"),
    string("ABLE"),
    string("PRI"),
    string("PRES"),
    string("ESE"),
    string("ANTI"),
    string("CALLY"),
    string("ATION"),
    string("EING"),
  };

class tpcw_worker : public bench_worker, public tpcw_worker_mixin {
public:
  tpcw_worker(unsigned int worker_id,
              unsigned long seed, ndb_wrapper *db,
              const map<string, ndb_ordered_index *> &open_tables,
              const map<string, vector<ndb_ordered_index *>> &partitions,
              spin_barrier *barrier_a, spin_barrier *barrier_b)
    : bench_worker(worker_id, true, seed, db,
                   open_tables, barrier_a, barrier_b),
      tpcw_worker_mixin(partitions)
  {
    if (verbose) {
      cerr << "tpcw: worker id " << worker_id << endl;
    }
    obj_key0.reserve(str_arena::MinStrReserveLength);
    obj_key1.reserve(str_arena::MinStrReserveLength);
    obj_v.reserve(str_arena::MinStrReserveLength);
    srand(worker_id);
    zeta_2_theta = zeta(2,gtheta);
    denom = zeta(NUM_ITEMS()-1,gtheta);
  }

    double zeta(uint64_t n, double theta) {
        double sum = 0;
        for (uint64_t i = 1; i <= n; i++)
                sum += pow(1.0 / i, theta);
        return sum;
    }

    uint64_t zipf(uint64_t n, double theta) {
        double alpha = 1 / (1 - theta);
        double zetan = denom;
        //double zetan = zeta(n,gtheta);
        double eta = (1 - pow(2.0 / n, 1 - theta)) /(1 - zeta_2_theta / zetan);
        //double eta = (1 - pow(2.0 / n, 1 - theta)) /(1 - zeta(2,gtheta) / zetan);
        double u = (double)(rand() % 10000000) / 10000000;
        double uz = u * zetan;
        if (uz < 1) return 1;
        if (uz < 1 + pow(0.5, theta)) return 2;
        return 1 + (uint64_t)(n * pow(eta*u -eta + 1, alpha));
    }


  txn_result txn_do_cart();

  static txn_result
  TxnDoCart(bench_worker *w)
  {
    return static_cast<tpcw_worker *>(w)->txn_do_cart();
  }

  txn_result txn_do_buy_confirm(); 

  static txn_result
  TxnDoBuyConfirm(bench_worker *w)
  {
    return static_cast<tpcw_worker *>(w)->txn_do_buy_confirm();
  }

  txn_result txn_do_buy_confirm2(); 
  static txn_result
  TxnDoBuyConfirm2(bench_worker *w)
  {
    return static_cast<tpcw_worker *>(w)->txn_do_buy_confirm2();
  }


  virtual workload_desc_vec
  get_workload() const
  {
    workload_desc_vec w;
    unsigned m = 0;
    for (size_t i = 0; i < ARRAY_NELEMS(g_txn_workload_mix); i++)
      m += g_txn_workload_mix[i];

    ALWAYS_ASSERT(m == 100);

    if(worker_id % coreid::static_cpus_online == 0){
        if(do_buy)
            w.push_back(workload_desc("DoBuyConfirm", 1, TxnDoBuyConfirm2));
        else{

            if (g_txn_workload_mix[0])
              w.push_back(workload_desc("DoCart", double(g_txn_workload_mix[0])/100.0, TxnDoCart));
            if (g_txn_workload_mix[1])
              w.push_back(workload_desc("DoBuyConfirm", double(g_txn_workload_mix[1])/100.0, TxnDoBuyConfirm));
        }
    }

    else{

        if (g_txn_workload_mix[0])
          w.push_back(workload_desc("DoCart", double(g_txn_workload_mix[0])/100.0, TxnDoCart));
        if (g_txn_workload_mix[1])
          w.push_back(workload_desc("DoBuyConfirm", double(g_txn_workload_mix[1])/100.0, TxnDoBuyConfirm));
    }


    return w;
  }

    int32_t get_next_id(int32_t id, int32_t index){
        int high = id>>9;
        int nid;
        high ++;
        nid = (high<<9) | index;
        return nid;
    }

    int32_t get_prev_id(int32_t id, int32_t index){
        int high = id>>9;
        int pid;
        high --;
        pid = (high<<9) | index;
        return pid;

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

  inline ALWAYS_INLINE string &
  str()
  {
    return *arena.next();
  }

private:

  // some scratch buffer space
  string obj_key0;
  string obj_key1;
  string obj_v;

  double denom;
  double zeta_2_theta;

};

/*
 * new txn template
 

tpcw_worker::txn_result
tpcw_worker::txn_do_cart()
{

    void *txn = db->new_txn(txn_flags, arena, txn_buf(), HINT_TPCW_DO_CART);

#ifdef USE_BCC
    bcc_txn_setup(txn);
#endif
#ifdef BCC_TXN_STATS
    bcc_txn_stats_setup(txn, "DO_CART");
#endif

    scoped_str_arena s_arena(arena);

    try{
        size_t ret = 0;

        measure_txn_counters(txn, "txn_do_cart");

        if (likely(db->commit_txn(txn,NULL))){
            return txn_result(true, ret);
        }else{
            db->abort_txn(txn);
            return txn_result(false, 0);
        }
    }catch(ndb_wrapper::ndb_abort_exception &ex) {
        db->abort_txn(txn);
    }

    return txn_result(false, 0);

}

*/



/*
 * Our implementation of do_cart is different from the spec.
 * We merge the create_empty_cart into the do_cart.
 *
 * 1. create a new shopping cart
 * 2. randomly add 1-5 items to the cart
 * 3. getCart and refresh Cart
 */

tpcw_worker::txn_result
tpcw_worker::txn_do_cart()
{
    int tindex = worker_id%coreid::static_cpus_online;

    void *txn = db->new_txn(txn_flags, arena, txn_buf(), HINT_TPCW_DO_CART);

#ifdef USE_BCC
    bcc_txn_setup(txn);
#endif
#ifdef BCC_TXN_STATS
    bcc_txn_stats_setup(txn, "DO_CART");
#endif

    scoped_str_arena s_arena(arena);

    try{

        string obj_buf;
        size_t ret = 0;

        // get the next shopping id

        next_id::key k_ni(tindex);
        next_id::value v_ni;

        tbl_next_id(tindex)->get(txn, Encode(k_ni), obj_v,-1);
        Decode(obj_v, v_ni);

        int shopping_id = v_ni.ni_shopping_id;

        v_ni.ni_shopping_id = get_next_id(shopping_id, tindex);
        tbl_next_id(tindex)->put(txn,Encode(str(), k_ni), Encode(str(), v_ni));

        // create a new shopping cart
        shopping_cart::key k_sc(shopping_id);
        shopping_cart::value v_sc(GetCurrentTimeMillis());

        tbl_shopping_cart(tindex)->put(txn,Encode(k_sc), Encode(obj_buf,v_sc));

        // add item to shopping cart_line
        shopping_cart_line::key k_scl;
        shopping_cart_line::value v_scl;

        k_scl.scl_sc_id = shopping_id;
        int i_num = 4; 

        for(int i=0;i<i_num;i++){
            //k_scl.scl_i_id = NonUniformRandom(r, 8191, 7911, 1, 100);
            k_scl.scl_i_id = RandomNumber(r, 1, NUM_ITEMS());
            v_scl.scl_qty = 1;

            // put: if tuple does not exist, it will insert the tuple
            tbl_shopping_cart_line(tindex)->put(txn,Encode(obj_key0,k_scl), Encode(obj_key1,v_scl));

        }

        // add random item to the cart if shopping cart is empty 
        // we add a "hot item" to the cart
        if((r.next()%100)<rperc)
            k_scl.scl_i_id = 1; 
        else
            k_scl.scl_i_id = RandomNumber(r, 1, NUM_ITEMS());
        v_scl.scl_qty = RandomNumber(r,1,5) ;

        // put: if tuple does not exist, it will insert the tuple
        tbl_shopping_cart_line(tindex)->put(txn,Encode(obj_key0,k_scl), Encode(obj_key1,v_scl));

        // refresh cart

        // reset cart time
        v_sc.sc_time = GetCurrentTimeMillis();
        tbl_shopping_cart(tindex)->put(txn,Encode(str(), k_sc), Encode(str(), v_sc));

        // get cart
 
        static_limit_callback<10> c(s_arena.get(), false); 
        const shopping_cart_line::key k_scl_0(shopping_id,1);
        const shopping_cart_line::key k_scl_1(shopping_id, NUM_ITEMS());

        tbl_shopping_cart_line(tindex)->scan(txn, Encode(obj_key0, k_scl_0), &Encode(obj_key1, k_scl_1), c, s_arena.get());

        ALWAYS_ASSERT(c.size()>0);

        for (size_t i = 0; i < c.size(); i++) {

            shopping_cart_line::key k_scl_temp;
            const shopping_cart_line::key *k_scl_p = Decode(*c.values[i].first, k_scl_temp);

            item::key k_i(k_scl_p->scl_i_id);
            ALWAYS_ASSERT(tbl_item(0)->get(txn, Encode(obj_key0, k_i), obj_v,-1));
        }

        measure_txn_counters(txn, "txn_do_cart");
        if (likely(db->commit_txn(txn,NULL))){
            return txn_result(true, ret);
        }else{
            db->abort_txn(txn);
            return txn_result(false, 0);
        }
    }catch(ndb_wrapper::ndb_abort_exception &ex) {
        db->abort_txn(txn);
    }

    return txn_result(false, 0);

}

/*
 * do_buy_confirm: submit a order based on the last shopping_id
 */


tpcw_worker::txn_result
tpcw_worker::txn_do_buy_confirm()
{

    int c_id = RandomNumber(r, 1, NUM_CUSTOMERS());
    int tindex = worker_id % coreid::static_cpus_online;

    void *txn = db->new_txn(txn_flags, arena, txn_buf(), HINT_TPCW_DO_BUY_CONFIRM);

#ifdef USE_BCC
    bcc_txn_setup(txn);
#endif
#ifdef BCC_TXN_STATS
    bcc_txn_stats_setup(txn, "DO_BUY_CONFIRM");
#endif

    scoped_str_arena s_arena(arena);

    try{
        size_t ret = 0;

        // get next id

        next_id::key k_ni(tindex);
        next_id::value v_ni;

        tbl_next_id(0)->get(txn, Encode(k_ni), obj_v,-1);
        Decode(obj_v, v_ni);

        int shopping_id = get_prev_id(v_ni.ni_shopping_id, tindex); 

        // get customer discount
        customer::key k_c(c_id);
        customer::value v_c;
        const customer::value * v_c_p;
        ALWAYS_ASSERT(tbl_customer(0)->get(txn, Encode(obj_key0, k_c), obj_v,-1));
        v_c_p = Decode(obj_v, v_c);
        float discount = v_c_p->c_discount;

        // get cart
        
        static_limit_callback<10> c(s_arena.get(), false); // never more than 5 items per shopping_cart 
        const shopping_cart_line::key k_scl_0(shopping_id,1);
        const shopping_cart_line::key k_scl_1(shopping_id, NUM_ITEMS());

        tbl_shopping_cart_line(0)->scan(txn, Encode(obj_key0, k_scl_0), &Encode(obj_key1, k_scl_1), c, s_arena.get());

        ALWAYS_ASSERT(c.size()>0);

        vector<int32_t> i_ids;
        vector<int32_t> i_qtys;

        for (size_t i = 0; i < c.size(); i++) {

            shopping_cart_line::key k_scl_temp;
            const shopping_cart_line::key *k_scl_p = Decode(*c.values[i].first, k_scl_temp);
            shopping_cart_line::value v_scl_temp;
            const shopping_cart_line::value *v_scl_p = Decode(*c.values[i].second, v_scl_temp);

            i_ids.push_back(k_scl_p->scl_i_id);
            i_qtys.push_back(v_scl_p->scl_qty);

            item::key k_i(k_scl_p->scl_i_id);
            ALWAYS_ASSERT(tbl_item(0)->get(txn, Encode(obj_key0, k_i), obj_v,-1));
        }

        // enter address
        int addr_id = v_ni.ni_addr_id; 

        v_ni.ni_addr_id = get_next_id(addr_id,tindex); 

        const address::key k_a(addr_id);

        address::value v_a;
        //const string addr_street1 = RandomStr(r, RandomNumber(r, 15, 40));
        v_a.addr_street1.assign("randomstreet");

        //const string addr_street2 = RandomStr(r, RandomNumber(r, 15, 40));
        v_a.addr_street2.assign("randomstreet");

        //const string addr_city = RandomStr(r, RandomNumber(r, 4,30));
        v_a.addr_city.assign("randomcity");

        //const string addr_state = RandomStr(r, RandomNumber(r, 2,30));
        v_a.addr_state.assign("randomstate");

        //const string addr_zip = RandomStr(r, RandomNumber(r, 5,10));
        v_a.addr_zip.assign("randomzip");

        v_a.addr_co_id = RandomNumber(r,1,92);

        tbl_address(0)->insert(txn, Encode(str(),k_a), Encode(str(), v_a));

        // enter order 
        int order_id = v_ni.ni_order_id;

        v_ni.ni_order_id = get_next_id(order_id,tindex); 

        tbl_next_id(0)->put(txn,Encode(str(), k_ni), Encode(str(), v_ni));


        orders::key k_o(order_id);
        orders::value v_o;

        v_o.o_c_id = c_id;
        v_o.o_date = GetCurrentTimeMillis();
        v_o.o_sub_total = (float) RandomNumber(r,1000,999999)/100.0 * discount;
        v_o.o_tax = v_o.o_sub_total * 0.0825;
        v_o.o_total = v_o.o_sub_total + v_o.o_tax + 3.00 + i_ids.size();
        v_o.o_ship_type.assign(ship_types[RandomNumber(r,0,num_ship_types-1)]);
        v_o.o_ship_date = GetCurrentTimeMillis();
        v_o.o_bill_addr_id = v_c_p->c_addr_id; 
        v_o.o_ship_addr_id = addr_id; 
        v_o.o_status.assign("Pending");

        tbl_orders(0)->insert(txn,Encode(str(),k_o), Encode(str(),v_o));

        int32_t id_len = i_ids.size();

        for(int32_t i=0;i<id_len;i++){
            int qty = i_qtys[i];

            order_line::key k_ol;
            order_line::value v_ol;
            
            k_ol.ol_id = i;
            k_ol.ol_o_id = order_id;
            v_ol.ol_i_id = i_ids.at(i); 
            v_ol.ol_qty = qty; 
            v_ol.ol_discount = discount; 
            //const string ol_comments = RandomStr(r, RandomNumber(r,20,100));
            //v_ol.ol_comments.assign(ol_comments);
            v_ol.ol_comments.assign("Random");

            tbl_order_line(0)->insert(txn, Encode(str(),k_ol), Encode(str(),v_ol));
            
            item::key k_i(i_ids.at(i));
            ALWAYS_ASSERT(tbl_item(0)->get(txn,Encode(obj_key0, k_i), obj_v, -1));
            item::value v_i_temp;
            Decode(obj_v, v_i_temp); 
            if(v_i_temp.i_stock - qty < 10){
                v_i_temp.i_stock = v_i_temp.i_stock - qty + 21;
            }else{
                v_i_temp.i_stock = v_i_temp.i_stock - qty;
            }
            tbl_item(0)->put(txn, Encode(str(), k_i), Encode(str(),v_i_temp));
        }

        // enter ccXacts

        const cc_xacts::key k_xact(order_id);
        cc_xacts::value v_xact;
        
        const string cx_type = credit_cards[RandomNumber(r,0, num_card_types-1)];
        v_xact.cx_type.assign(cx_type);

        //const string cx_num = RandomStr(r,16);
        v_xact.cx_num.assign("randomnumber");

        //const string cx_name = RandomStr(r, RandomNumber(r,14,30));
        v_xact.cx_name.assign("randomename");

        v_xact.cx_expiry = GetCurrentTimeMillis();
        //const string cx_auth_id = RandomStr(r, 15);
        v_xact.cx_auth_id.assign("randomid");
        v_xact.cx_co_id = RandomNumber(r,1,92);

        tbl_cc_xacts(0)->insert(txn,Encode(str(),k_xact), Encode(str(),v_xact));

        // clear cart

        for(int32_t i=0;i<id_len;i++){
            shopping_cart_line::key k_sc(shopping_id,i_ids.at(i));
            //tbl_shopping_cart_line(0)->remove(txn, Encode(str(), k_sc));
        }

        measure_txn_counters(txn, "txn_do_buy_confirm");

        if (likely(db->commit_txn(txn,NULL))){
            return txn_result(true, ret);
        }else{
            db->abort_txn(txn);
            return txn_result(false, 0);
        }
    }catch(ndb_wrapper::ndb_abort_exception &ex) {
        db->abort_txn(txn);
    }

    return txn_result(false, 0);

}

/*
 * do_buy_confirm2: simply update the item table.
 */


tpcw_worker::txn_result
tpcw_worker::txn_do_buy_confirm2()
{

    void *txn = db->new_txn(txn_flags, arena, txn_buf(), HINT_TPCW_DO_BUY_CONFIRM);

#ifdef USE_BCC
    bcc_txn_setup(txn);
#endif
#ifdef BCC_TXN_STATS
    bcc_txn_stats_setup(txn, "DO_BUY_CONFIRM2");
#endif

    scoped_str_arena s_arena(arena);

    try{
        size_t ret = 0;

        //int i_num = RandomNumber(r,1,5);
        int i_num = 1; 

        for(int32_t i=0;i<i_num;i++){
            int qty = RandomNumber(r,1,5);
            //int i_id = zipf(NUM_ITEMS()-1, gtheta); 
            int i_id = i+1; 

            item::key k_i(i_id);
            ALWAYS_ASSERT(tbl_item(0)->get(txn,Encode(obj_key0, k_i), obj_v, -1));
            item::value v_i_temp;
            Decode(obj_v, v_i_temp);
            if(v_i_temp.i_stock - qty < 10){
                v_i_temp.i_stock = v_i_temp.i_stock - qty + 21;
            }else{
                v_i_temp.i_stock = v_i_temp.i_stock - qty;
            }
            tbl_item(0)->put(txn, Encode(str(), k_i), Encode(str(),v_i_temp));
        }
        measure_txn_counters(txn, "txn_do_buy_confirm2");

        if (likely(db->commit_txn(txn,NULL))){
            return txn_result(true, ret);
        }else{
            db->abort_txn(txn);
            return txn_result(false, 0);
        }
    }catch(ndb_wrapper::ndb_abort_exception &ex) {
        db->abort_txn(txn);
    }

    return txn_result(false, 0);

}



class tpcw_next_id_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_next_id_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {
        
    zeta_2_theta = zeta(2,gtheta);
    denom = zeta(NUM_ITEMS()-1,gtheta);
  }


    double zeta(uint64_t n, double theta) {
        double sum = 0;
        for (uint64_t i = 1; i <= n; i++)
                sum += pow(1.0 / i, theta);
        return sum;
    }

    uint64_t zipf(uint64_t n, double theta) {
        double alpha = 1 / (1 - theta);
        double zetan = denom;
        //double zetan = zeta(n,gtheta);
        double eta = (1 - pow(2.0 / n, 1 - theta)) /(1 - zeta_2_theta / zetan);
        //double eta = (1 - pow(2.0 / n, 1 - theta)) /(1 - zeta(2,gtheta) / zetan);
        double u = (double)(rand() % 10000000) / 10000000;
        double uz = u * zetan;
        if (uz < 1) return 1;
        if (uz < 1 + pow(0.5, theta)) return 2;
        return 1 + (uint64_t)(n * pow(eta*u -eta + 1, alpha));
    }



protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 0; i <coreid::static_cpus_online; i++) {

        int ni_id;

        for(int j=0;j<NUM_SCL();j++){
            ni_id = (j<<9) | i;
            //int i_num = RandomNumber(r,1,5);
            int i_num = 1; 

            for(int k=0;k<i_num;k++){
                int i_qty = RandomNumber(r,1,5);
                int i_id = zipf(NUM_ITEMS()-1, gtheta); 

                const shopping_cart_line::key k_scl (ni_id,i_id);
                const shopping_cart_line::value v_scl (i_qty);

                tbl_shopping_cart_line(i)->insert(txn, Encode(k_scl), Encode(obj_buf,v_scl));
            }
        }

        ni_id = ((NUM_SCL())<<9)|i;
        const next_id::key k_ni(i);
        const next_id::value v_ni(ni_id,ni_id,ni_id);
        tbl_next_id(i)->insert(txn, Encode(k_ni), Encode(obj_buf, v_ni));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading next_id" << endl;
      cerr << "[INFO]   * average next_id record length: "
           << (double(total_sz)/double(NUM_IDS())) << " bytes" << endl;
    }
  }
private:

  double denom;
  double zeta_2_theta;

};



class tpcw_orderxact_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_orderxact_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_ORDERS(); i++) {

        const orders::key k(i);
        orders::value v;

        int num_items = RandomNumber(r,1,5); 
        v.o_c_id = RandomNumber(r,1,NUM_CUSTOMERS());
        v.o_date = GetCurrentTimeMillis(); 
        v.o_sub_total = (float) RandomNumber(r,1000,999999)/100.0;
        v.o_tax = v.o_sub_total * 0.0825;
        v.o_total = v.o_sub_total + v.o_tax + 3.00 + num_items;
        v.o_ship_type.assign(ship_types[RandomNumber(r,0,num_ship_types-1)]);
        v.o_ship_date = GetCurrentTimeMillis(); 
        v.o_bill_addr_id = RandomNumber(r,1, 2*NUM_CUSTOMERS());
        v.o_ship_addr_id = RandomNumber(r,1, 2*NUM_CUSTOMERS());
        v.o_status.assign(status_types[RandomNumber(r,0,num_status_types-1)]);

        tbl_orders(0)->insert(txn,Encode(k), Encode(obj_buf,v));

        for(int j = 1; j <= num_items; j++){
            order_line::key kk;
            order_line::value vv;
            
            kk.ol_id = j;
            kk.ol_o_id = i;
            vv.ol_i_id = RandomNumber(r,1, NUM_ITEMS());
            vv.ol_qty = RandomNumber(r,1, 300);
            vv.ol_discount = RandomNumber(r,0,30)/100.0;
            const string ol_comments = RandomStr(r, RandomNumber(r,20,100));
            vv.ol_comments.assign(ol_comments);

            tbl_order_line(0)->insert(txn, Encode(kk), Encode(obj_buf,vv));
         }

        const cc_xacts::key kkk(i);
        cc_xacts::value vvv;
        
        const string cx_type = credit_cards[RandomNumber(r,0, num_card_types-1)];
        vvv.cx_type.assign(cx_type);

        const string cx_num = RandomStr(r,16);
        vvv.cx_num.assign(cx_num);

        const string cx_name = RandomStr(r, RandomNumber(r,14,30));
        vvv.cx_name.assign(cx_name);

        vvv.cx_expiry = GetCurrentTimeMillis();
        const string cx_auth_id = RandomStr(r, 15);
        vvv.cx_auth_id.assign(cx_auth_id);
        vvv.cx_co_id = RandomNumber(r,1,92);

        tbl_cc_xacts(0)->insert(txn,Encode(kkk), Encode(obj_buf,vvv));

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_author(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading orders and cx_xacts" << endl;
    }
  }
};




class tpcw_country_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_country_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string countries[] = {"United States","United Kingdom","Canada",
                              "Germany","France","Japan","Netherlands",
                              "Italy","Switzerland", "Australia","Algeria",
                              "Argentina","Armenia","Austria","Azerbaijan",
                              "Bahamas","Bahrain","Bangla Desh","Barbados",
                              "Belarus","Belgium","Bermuda", "Bolivia",
                              "Botswana","Brazil","Bulgaria","Cayman Islands",
                              "Chad","Chile", "China","Christmas Island",
                              "Colombia","Croatia","Cuba","Cyprus",
                              "Czech Republic","Denmark","Dominican Republic",
                              "Eastern Caribbean","Ecuador", "Egypt",
                              "El Salvador","Estonia","Ethiopia",
                              "Falkland Island","Faroe Island", "Fiji",
                              "Finland","Gabon","Gibraltar","Greece","Guam",
                              "Hong Kong","Hungary", "Iceland","India",
                              "Indonesia","Iran","Iraq","Ireland","Israel",
                              "Jamaica", "Jordan","Kazakhstan","Kuwait",
                              "Lebanon","Luxembourg","Malaysia","Mexico",
                              "Mauritius", "New Zealand","Norway","Pakistan",
                              "Philippines","Poland","Portugal","Romania",
                              "Russia","Saudi Arabia","Singapore","Slovakia",
                              "South Africa","South Korea", "Spain","Sudan",
                              "Sweden","Taiwan","Thailand","Trinidad",
                              "Turkey","Venezuela", "Zambia"};

        float exchanges[] = { 1, .625461, 1.46712, 1.86125, 6.24238, 121.907,
                               2.09715, 1842.64, 1.51645, 1.54208, 65.3851,
                               0.998, 540.92, 13.0949, 3977, 1, .3757,
                               48.65, 2, 248000, 38.3892, 1, 5.74, 4.7304,
                               1.71, 1846, .8282, 627.1999, 494.2, 8.278,
                               1.5391, 1677, 7.3044, 23, .543, 36.0127,
                               7.0707, 15.8, 2.7, 9600, 3.33771, 8.7,
                               14.9912, 7.7, .6255, 7.124, 1.9724, 5.65822,
                               627.1999, .6255, 309.214, 1, 7.75473, 237.23,
                               74.147, 42.75, 8100, 3000, .3083, .749481,
                               4.12, 37.4, 0.708, 150, .3062, 1502, 38.3892,
                               3.8, 9.6287, 25.245, 1.87539, 7.83101,
                               52, 37.8501, 3.9525, 190.788, 15180.2,
                               24.43, 3.7501, 1.72929, 43.9642, 6.25845,
                               1190.15, 158.34, 5.282, 8.54477, 32.77, 37.1414,
                               6.1764, 401500, 596, 2447.7 };

    string currencies[] = { "Dollars","Pounds","Dollars","Deutsche Marks",
                                "Francs","Yen","Guilders","Lira","Francs",
                                "Dollars","Dinars","Pesos", "Dram",
                                "Schillings","Manat","Dollars","Dinar","Taka",
                                "Dollars","Rouble","Francs","Dollars",
                                "Boliviano", "Pula", "Real", "Lev","Dollars",
                                "Franc","Pesos","Yuan Renmimbi","Dollars",
                                "Pesos","Kuna","Pesos","Pounds","Koruna",
                                "Kroner","Pesos","Dollars","Sucre","Pounds",
                                "Colon","Kroon","Birr","Pound","Krone",
                                "Dollars","Markka","Franc","Pound","Drachmas",
                                "Dollars","Dollars","Forint","Krona","Rupees",
                                "Rupiah","Rial","Dinar","Punt","Shekels",
                                "Dollars","Dinar","Tenge","Dinar","Pounds",
                                "Francs","Ringgit","Pesos","Rupees","Dollars",
                                "Kroner","Rupees","Pesos","Zloty","Escudo",
                                "Leu","Rubles","Riyal","Dollars","Koruna",
                                "Rand","Won","Pesetas","Dinar","Krona",
                                "Dollars","Baht","Dollars","Lira","Bolivar",
                                "Kwacha"};

    uint NUM_COUNTRIES = 92;

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_COUNTRIES; i++) {

        const country::key k(i);

        country::value v;

        v.co_name.assign(countries[i-1]);
        v.co_exchange = exchanges[i-1];
        v.co_currency.assign(currencies[i-1]);

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_country(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading country" << endl;
      cerr << "[INFO]   * average item record length: "
           << (double(total_sz)/double(NUM_COUNTRIES)) << " bytes" << endl;
    }
  }
};



class tpcw_author_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_author_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_AUTHORS(); i++) {

        const author::key k(i);

        author::value v;

        const string a_fname = RandomStr(r, RandomNumber(r,3,20));
        v.a_fname.assign(a_fname);

        const string a_mname = RandomStr(r, RandomNumber(r,1,20));
        v.a_mname.assign(a_mname);

        const string a_lname = RandomStr(r, RandomNumber(r,1,20));
        v.a_lname.assign(a_lname);

        v.a_dob = GetCurrentTimeMillis();

        const string a_bio = RandomStr(r, RandomNumber(r,125,500));
        v.a_bio.assign(a_bio);

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_author(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading author" << endl;
      cerr << "[INFO]   * average item record length: "
           << (double(total_sz)/double(NUM_AUTHORS())) << " bytes" << endl;
    }
  }
};




class tpcw_address_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_address_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_ADDRESSES(); i++) {

        const address::key k(i);

        address::value v;
        const string addr_street1 = RandomStr(r, RandomNumber(r, 15, 40));
        v.addr_street1.assign(addr_street1);

        const string addr_street2 = RandomStr(r, RandomNumber(r, 15, 40));
        v.addr_street2.assign(addr_street2);

        const string addr_city = RandomStr(r, RandomNumber(r, 4,30));
        v.addr_city.assign(addr_city);

        const string addr_state = RandomStr(r, RandomNumber(r, 2,30));
        v.addr_state.assign(addr_state);

        const string addr_zip = RandomStr(r, RandomNumber(r, 5,10));
        v.addr_zip.assign(addr_zip);

        v.addr_co_id = RandomNumber(r,1,92);

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_address(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading address" << endl;
      cerr << "[INFO]   * average item record length: "
           << (double(total_sz)/double(NUM_ADDRESSES())) << " bytes" << endl;
    }
  }
};



class tpcw_customer_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_customer_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_CUSTOMERS(); i++) {

        const customer::key k(i);

        customer::value v;
        const string c_uname = RandomStr(r, RandomNumber(r, 10, 20));
        v.c_uname.assign(c_uname);

        const string c_passwd = RandomStr(r, RandomNumber(r, 10, 20));
        v.c_passwd.assign(c_passwd);

        const string c_lname = RandomStr(r, RandomNumber(r, 8, 15));
        v.c_lname.assign(c_lname);

        const string c_fname = RandomStr(r, RandomNumber(r, 8, 15));
        v.c_fname.assign(c_fname);

        v.c_addr_id = RandomNumber(r,1, 2* NUM_CUSTOMERS());

        const string c_phone = RandomStr(r, RandomNumber(r, 9, 16));
        v.c_phone.assign(c_phone);

        const string c_email = RandomStr(r, RandomNumber(r, 14, 30));
        v.c_email.assign(c_email);

        v.c_since = GetCurrentTimeMillis();

        v.c_last_login = GetCurrentTimeMillis();

        v.c_login =  GetCurrentTimeMillis();

        v.c_expiration = GetCurrentTimeMillis();

        v.c_discount = (float) RandomNumber(r,0,50)/100.0;

        v.c_balance = 0;

        v.c_ytd_pmt = (float) RandomNumber(r,0,99999)/100.0;

        v.c_birthdate = GetCurrentTimeMillis();

        const string c_data = RandomStr(r, RandomNumber(r,100,500));
        v.c_data.assign(c_data);

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_customer(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading customer" << endl;
      cerr << "[INFO]   * average item record length: "
           << (double(total_sz)/double(NUM_CUSTOMERS())) << " bytes" << endl;
    }
  }
};



class tpcw_item_loader : public bench_loader, public tpcw_worker_mixin {
public:
  tpcw_item_loader(unsigned long seed,
                   ndb_wrapper *db,
                   const map<string, ndb_ordered_index *> &open_tables,
                   const map<string, vector<ndb_ordered_index *>> &partitions)
    : bench_loader(seed, db, open_tables),
      tpcw_worker_mixin(partitions)
  {}

protected:
  virtual void
  load()
  {
#ifdef USE_BCC
    uint64_t latest_tid = 0, tid;
#endif

    string obj_buf;

    const ssize_t bsize = db->txn_max_batch_size();
    void *txn = db->new_txn(txn_flags, arena, txn_buf());

    uint64_t total_sz = 0;
    try {
      for (uint i = 1; i <= NUM_ITEMS(); i++) {
        // items don't "belong" to a certain warehouse, so no pinning
        const item::key k(i);

        item::value v;
        const string i_title = RandomStr(r, RandomNumber(r, 14, 60));
        v.i_title.assign(i_title);

        if(i<= NUM_ITEMS()/4)
          v.i_a_id = i;
        else
          v.i_a_id = RandomNumber(r, 1, NUM_ITEMS()/4);

        v.i_pub_date =  GetCurrentTimeMillis();
        const string i_publisher = RandomStr(r, RandomNumber(r, 14, 60));
        v.i_publisher.assign(i_publisher);

        const string i_subject = RandomStr(r, RandomNumber(r, 14, 60));
        v.i_subject.assign(i_subject);

        const string i_desc = RandomStr(r, RandomNumber(r, 100, 500));
        v.i_desc.assign(i_desc);

        v.i_related1 = RandomNumber(r,1, NUM_ITEMS());

        do{
            v.i_related2 = RandomNumber(r, 1, NUM_ITEMS());
        }while(v.i_related1 == v.i_related2);

        do{
            v.i_related3 = RandomNumber(r,1, NUM_ITEMS());
        }while(v.i_related1 == v.i_related3 || v.i_related2 == v.i_related3);

        do{
            v.i_related4 = RandomNumber(r,1, NUM_ITEMS());
        }while(v.i_related1 == v.i_related4 || v.i_related2 == v.i_related4 || v.i_related3 == v.i_related4);

        do{
            v.i_related5 = RandomNumber(r,1, NUM_ITEMS());
        }while(v.i_related1 == v.i_related5 || v.i_related2 == v.i_related5 || v.i_related3 == v.i_related5 || v.i_related4 == v.i_related5);

        const string i_thumbnail = RandomStr(r, RandomNumber(r,14,60));
        v.i_thumbnail.assign(i_thumbnail);

        const string i_image = RandomStr(r, RandomNumber(r,14,60));
        v.i_image.assign(i_image);

        v.i_srp = (float) RandomNumber(r, 100,99999)/100.0;

        v.i_cost = v.i_srp - (RandomNumber(r,0,50)/100.0)*v.i_srp;

        v.i_avail = GetCurrentTimeMillis();

        v.i_stock = RandomNumber(r,10,30);
        v.i_isbn = RandomStr(r, 12);
        v.i_page = RandomNumber(r,20,9999);

        const string i_backing = RandomStr(r, RandomNumber(r, 1,15));
        v.i_backing.assign(i_backing);

        const string i_dimensions = RandomStr(r,RandomNumber(r, 1,25));
        v.i_dimensions.assign(i_dimensions);

        const size_t sz = Size(v);
        total_sz += sz;
        tbl_item(0)->insert(txn, Encode(k), Encode(obj_buf, v));

        if (bsize != -1 && !(i % bsize)) {
#ifdef USE_BCC
          ALWAYS_ASSERT(db->commit_txn(txn, &tid));
          if(tid > latest_tid)
            latest_tid = tid;
#else
          ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
          txn = db->new_txn(txn_flags, arena, txn_buf());
          arena.reset();
        }
      }
#ifdef USE_BCC
      ALWAYS_ASSERT(db->commit_txn(txn, &tid));
      if(tid > latest_tid)
        latest_tid = tid;
#else
      ALWAYS_ASSERT(db->commit_txn(txn,NULL));
#endif
    } catch (ndb_wrapper::ndb_abort_exception &ex) {
      // shouldn't abort on loading!
      ALWAYS_ASSERT(false);
    }
#ifdef USE_BCC
    loader_tids_[loader_id_] = latest_tid;
#endif
    if (verbose) {
      cerr << "[INFO] finished loading item" << endl;
      cerr << "[INFO]   * average item record length: "
           << (double(total_sz)/double(NUM_ITEMS())) << " bytes" << endl;
    }
  }
};


template <typename T>
static vector<T>
unique_filter(const vector<T> &v)
{
  set<T> seen;
  vector<T> ret;
  for (auto &e : v)
    if (!seen.count(e)) {
      ret.emplace_back(e);
      seen.insert(e);
    }
  return ret;
}

class tpcw_bench_runner : public bench_runner {
private:

  static vector<ndb_ordered_index *>
  OpenTablesForTablespace(ndb_wrapper *db, const char *name, size_t expected_size)
  {
    const string s_name(name);
    int n = coreid::num_cpus_online();
    vector<ndb_ordered_index *> ret(n);
    for(int i=0;i<n;i++)
        ret[i] = db->open_index(s_name + "_" + to_string(i), expected_size, false);
    return ret;
  }

public:
  tpcw_bench_runner(ndb_wrapper *db, size_t nthreads)
    : bench_runner(db, nthreads)
  {

#define OPEN_TABLESPACE_X(x) \
    partitions[#x] = OpenTablesForTablespace(db, #x, sizeof(x));

    TPCW_TABLE_LIST(OPEN_TABLESPACE_X);

#undef OPEN_TABLESPACE_X

    for (auto &t : partitions) {
      auto v = unique_filter(t.second);
      for (size_t i = 0; i < v.size(); i++)
        open_tables[t.first + "_" + to_string(i)] = v[i];
    }
  }

protected:
  virtual vector<bench_loader *>
  make_loaders()
  {
    vector<bench_loader *> ret;
    ret.push_back(new tpcw_item_loader(235443, db, open_tables, partitions));
    ret.push_back(new tpcw_next_id_loader(1, db, open_tables, partitions));
    ret.push_back(new tpcw_customer_loader(89785943, db, open_tables, partitions));
    return ret;
  }

  virtual vector<bench_worker *>
  make_workers()
  {
    const unsigned alignment = coreid::num_cpus_online();
    const int blockstart =
      coreid::allocate_contiguous_aligned_block(nthreads, alignment);
    ALWAYS_ASSERT(blockstart >= 0);
    ALWAYS_ASSERT((blockstart % alignment) == 0);
    fast_random r(23984543);
    vector<bench_worker *> ret;
      for (size_t i = 0; i < nthreads; i++)
        ret.push_back(
          new tpcw_worker(
            blockstart + i,
            r.next(), db, open_tables, partitions,
            &barrier_a, &barrier_b));

    return ret;
  }

private:
  map<string, vector<ndb_ordered_index *>> partitions;
};

void
bcctest(ndb_wrapper *db, size_t nthreads, int argc, char **argv)
{

  // parse options
  optind = 1;
  while (1) {
    static struct option long_options[] =
    {
      {"workload-mix"                         , required_argument , 0                                     , 'w'} ,
      {"read-percent"                         , required_argument , 0                                     , 'r'} ,
      {"gtheta"                         , required_argument , 0                                     , 'g'} ,
      {"--disable-do-buy"                         , required_argument , 0                                     , 'd'} ,
      {0, 0, 0, 0}
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "r:w:g:d", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {

    case 'd':
        do_buy = 0; 
        break;

    case 'g':
        gtheta = atoi(optarg)/100.0;
        break;

    case 'r':
        rperc = atoi(optarg);
        break;

    case 'w':
      {
        const vector<string> toks = split(optarg, ',');
        ALWAYS_ASSERT(toks.size() == ARRAY_NELEMS(g_txn_workload_mix));
        unsigned s = 0;
        for (size_t i = 0; i < toks.size(); i++) {
          unsigned p = strtoul(toks[i].c_str(), nullptr, 10);
          ALWAYS_ASSERT(p >= 0 && p <= 100);
          s += p;
          g_txn_workload_mix[i] = p;
        }
        ALWAYS_ASSERT(s == 100);
      }
      break;

    case '?':
      /* getopt_long already printed an error message. */
      exit(1);

    default:
      abort();
    }
  }


  if (verbose) {
    cerr << "tpcw settings:" << endl;
    cerr << "  workload_mix                 : " <<
      format_list(g_txn_workload_mix,
                  g_txn_workload_mix + ARRAY_NELEMS(g_txn_workload_mix)) << endl;
  }

  tpcw_bench_runner r(db, nthreads);
  r.run();
  exit(0); 
}
