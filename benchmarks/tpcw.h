#ifndef _NDB_BENCH_TPCW_H_
#define _NDB_BENCH_TPCW_H_

#include "../record/encoder.h"
#include "../record/inline_str.h"
#include "../macros.h"

#define ADDRESS_KEY_FIELDS(x,y) \
  x(int32_t, addr_id)

#define ADDRESS_VALUE_FILEDS(x,y)  \
  x(inline_str_8<40>,addr_street1)  \
  y(inline_str_8<40>,addr_street2)  \
  y(inline_str_8<30>,addr_city)  \
  y(inline_str_8<20>,addr_state)  \
  y(inline_str_8<10>,addr_zip)  \
  y(int32_t,addr_co_id)
DO_STRUCT(address, ADDRESS_KEY_FIELDS, ADDRESS_VALUE_FILEDS)

#define AUTHOR_KEY_FIELDS(x,y)  \
  x(int32_t, a_id)

#define AUTHOR_VALUE_FIELDS(x,y) \
  x(inline_str_8<20>, a_fname)  \
  y(inline_str_8<20>, a_lname)  \
  y(inline_str_8<20>, a_mname)  \
  y(int32_t, a_dob)  \
  y(inline_str_8<500>, a_bio) 
DO_STRUCT(author, AUTHOR_KEY_FIELDS, AUTHOR_VALUE_FIELDS)

#define CC_XACTS_KEY_FIELDS(x,y)  \
  x(int32_t, cx_o_id)
#define CC_XACTS_VALUE_FIELDS(x,y)  \
  x(inline_str_8<10>, cx_type)  \
  y(inline_str_8<20>, cx_num)  \
  y(inline_str_8<30>, cx_name)  \
  y(int32_t, cx_expiry)  \
  y(inline_str_8<15>, cx_auth_id)  \
  y(float, cx_xact_amt)  \
  y(int32_t, cx_xact_date)  \
  y(int32_t, cx_co_id) 
DO_STRUCT(cc_xacts, CC_XACTS_KEY_FIELDS,CC_XACTS_VALUE_FIELDS)

#define COUNTRY_KEY_FIELDS(x,y) \
  x(int32_t, co_id)
#define COUNTRY_VALUE_FIELDS(x,y)  \
  x(inline_str_8<50>, co_name)  \
  y(float,co_exchange)  \
  y(inline_str_8<18>, co_currency)
DO_STRUCT(country, COUNTRY_KEY_FIELDS, COUNTRY_VALUE_FIELDS)
  

#define CUSTOMER_KEY_FIELDS(x,y) \
  x(int32_t, c_id)

#define CUSTOMER_VALUE_FIELDS(x,y) \
  x(inline_str_8<20>,c_uname) \
  y(inline_str_8<20>,c_passwd) \
  y(inline_str_8<17>,c_fname) \
  y(inline_str_8<17>,c_lname) \
  y(int32_t, c_addr_id) \
  y(inline_str_8<18>, c_phone)  \
  y(inline_str_8<50>, c_email)  \
  y(int32_t, c_since) \
  y(int32_t, c_last_login)  \
  y(int32_t, c_login) \
  y(int32_t, c_expiration)  \
  y(float, c_discount)  \
  y(float, c_balance) \
  y(float, c_ytd_pmt) \
  y(int32_t, c_birthdate) \
  y(inline_str_8<510>, c_data)
DO_STRUCT(customer,CUSTOMER_KEY_FIELDS,CUSTOMER_VALUE_FIELDS)

#define ITEM_KEY_FIELDS(x, y) \
  x(int32_t,i_id)
#define ITEM_VALUE_FIELDS(x, y) \
  x(inline_str_8<60>,i_title) \
  y(int32_t, i_a_id)  \
  y(int32_t, i_pub_date)  \
  y(inline_str_8<60>, i_publisher)  \
  y(inline_str_8<60>, i_subject)  \
  y(inline_str_8<500>, i_desc)  \
  y(int32_t, i_related1)  \
  y(int32_t, i_related2)  \
  y(int32_t, i_related3)  \
  y(int32_t, i_related4)  \
  y(int32_t, i_related5)  \
  y(inline_str_8<40>, i_thumbnail)  \
  y(inline_str_8<40>, i_image)  \
  y(float, i_srp)  \
  y(float,i_cost) \
  y(int32_t,i_avail)  \
  y(int32_t, i_stock) \
  y(inline_str_8<12>, i_isbn) \
  y(int32_t, i_page)  \
  y(inline_str_8<15>, i_backing)  \
  y(inline_str_8<25>, i_dimensions)
DO_STRUCT(item, ITEM_KEY_FIELDS, ITEM_VALUE_FIELDS)

#define ORDER_LINE_KEY_FIELDS(x, y) \
  x(int32_t,ol_id) \
  y(int32_t,ol_o_id) 

#define ORDER_LINE_VALUE_FIELDS(x, y) \
  x(int32_t,ol_i_id) \
  y(int32_t,ol_qty) \
  y(float,ol_discount) \
  y(inline_str_8<110>, ol_comments)
DO_STRUCT(order_line, ORDER_LINE_KEY_FIELDS, ORDER_LINE_VALUE_FIELDS)

#define ORDERS_KEY_FIELDS(x,y)  \
  x(int32_t, o_id)
#define ORDERS_VALUE_FIELDS(x,y)  \
  x(int32_t, o_c_id)  \
  y(int32_t, o_date)  \
  y(float, o_sub_total) \
  y(float, o_tax) \
  y(float, o_total) \
  y(inline_str_8<20>, o_ship_type)  \
  y(int32_t, o_ship_date) \
  y(int32_t, o_bill_addr_id)  \
  y(int32_t, o_ship_addr_id)  \
  y(inline_str_8<20>, o_status)
DO_STRUCT(orders, ORDERS_KEY_FIELDS, ORDERS_VALUE_FIELDS)

#define SHOPPING_CART_KEY_FIELDS(x,y) \
  x(int32_t, sc_id) 
#define SHOPPING_CART_VALUE_FIELDS(x,y) \
  x(int32_t, sc_time)
DO_STRUCT(shopping_cart,SHOPPING_CART_KEY_FIELDS,SHOPPING_CART_VALUE_FIELDS)

#define SHOPPING_CART_LINE_KEY_FIELDS(x,y)  \
  x(int32_t, scl_sc_id) \
  y(int32_t, scl_i_id)
#define SHOPPING_CART_LINE_VALUE_FIELDS(x,y)  \
  x(int32_t, scl_qty)
DO_STRUCT(shopping_cart_line,SHOPPING_CART_LINE_KEY_FIELDS, SHOPPING_CART_LINE_VALUE_FIELDS)

#define NEXT_ID_KEY_FIELDS(x,y) \
  x(int32_t,ni_id)
#define NEXT_ID_VALUE_FIELDS(x,y) \
  x(int32_t,ni_shopping_id) \
  y(int32_t,ni_order_id)  \
  y(int32_t,ni_addr_id)
DO_STRUCT(next_id,NEXT_ID_KEY_FIELDS, NEXT_ID_VALUE_FIELDS)

#endif
