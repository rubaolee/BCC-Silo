#ifndef __CC_SCHEMA__
#define __CC_SCHEMA__

#include <vector>
#include "record/encoder.h"

#define SCAN_KEY(x,y)   \
        x(uint32_t,t1)  

#define SCAN_VALUE(x,y) \
        x(uint32_t,t2)

DO_STRUCT(testscan, SCAN_KEY, SCAN_VALUE)


/*
 * Table index.
 */

enum{
    TESTSCAN = 0,

    CUSTOMER,
    CUSTOMER_NAME_IDX,
    DISTRICT,
    HISTORY,
    ITEM,
    NEW_ORDER,
    OORDER,
    OORDER_C_ID_IDX,
    ORDER_LINE,
    STOCK,
    STOCK_DATA,
    WAREHOUSE
};


class my_scan_callback{
//class my_scan_callback : public abstract_ordered_index::new_scan_callback{

public:
    my_scan_callback(uint32_t min_v, uint32_t max_v):min(min_v),max(max_v){}
    std::vector<testscan::value> res;

    inline bool invoke(const std::string &key,   const uint8_t *value){

        testscan::value vv;
        const testscan::value *cvv;

        cvv = Decode(value,vv);
        if(cvv->t2 >=min && cvv->t2<= max )
            res.push_back(*cvv);

        return true;
    }

private:
    uint32_t min;
    uint32_t max;

};


#endif
