#ifndef __SH_HT_H__
#define __SH_HT_H__

#include <cstdint>
#ifdef BCC_SEARCH_TUPLE_SSE
#include <smmintrin.h>
#endif

#include "macros.h"
#include "sh_macros.h"

class TupleSearchInterface {
public:
  virtual ~TupleSearchInterface() {}

  // Checks whether tuple exists in the tuple set.
  virtual bool SearchTuple(const void **tuples, int ntuples/*const void *tuple*/) const = 0;
#ifdef BCC_SHHT_STATS
  // Returns the fullness of the tuple buckets.
  virtual double Fullness() const = 0;
#endif
};

template <size_t SetSize, size_t BucketsSize>
class HashedTupleSet : public TupleSearchInterface {
protected:
  struct TupleSetElem {
    const void *tuple_;
    short hash_next_;
  };

public:
  HashedTupleSet() : n_(0)
#ifdef BCC_SHHT_STATS
    , n_buckets_occupied_(0)
#endif
  {
    for (size_t i = 0; i < BucketsSize; ++i)
      hash_first_[i] = -1;
  }

  virtual bool SearchTuple(const void **tuples, int ntuples/*const void *tuple*/) const OVERRIDE {
    for (int it = 0; it < ntuples; ++it) {
      const uint64_t key = HashKey(tuples[it]);
      for (short i = hash_first_[key]; i != -1; i = tuple_set_[i].hash_next_) {
        if (tuple_set_[i].tuple_ == tuples[it])
          return true;
      }
    }
    return false;
  }

#ifdef BCC_SHHT_STATS
  virtual double Fullness() const OVERRIDE {
    return static_cast<double>(n_buckets_occupied_) / BucketsSize;
  }
#endif

  inline void AddTuple(const void *tuple) {
#ifdef BCC_CHECK
    ALWAYS_ASSERT(n_ < SetSize);
#endif
    const uint64_t key = HashKey(tuple);
    tuple_set_[n_].tuple_ = tuple;
    tuple_set_[n_].hash_next_ = hash_first_[key];
#ifdef BCC_SHHT_STATS
    if (hash_first_[key] == -1)
      ++n_buckets_occupied_;
#endif
    hash_first_[key] = n_;
    ++n_;
  }

protected:
  // Hard-coded hash function for various BucketsSize. The selection
  // of bits is based on offline bit-significance analysis.
  #define HASHKEY5(val) ( \
    ((val & (1UL << 11)) >> 11) | \
    ((val & (1UL << 8)) >> 7) | \
    ((val & (1UL << 14)) >> 12) | \
    ((val & (1UL << 10)) >> 7) | \
    ((val & (1UL << 9)) >> 5))

  #define HASHKEY6(val) ( \
    ((val & (1UL << 11)) >> 11) | \
    ((val & (1UL << 8)) >> 7) | \
    ((val & (1UL << 14)) >> 12) | \
    ((val & (1UL << 10)) >> 7) | \
    ((val & (1UL << 9)) >> 5) | \
    ((val & (1UL << 27)) >> 22))

  #define HASHKEY7(val) ( \
    ((val & (1UL << 11)) >> 11) | \
    ((val & (1UL << 8)) >> 7) | \
    ((val & (1UL << 14)) >> 12) | \
    ((val & (1UL << 10)) >> 7) | \
    ((val & (1UL << 9)) >> 5) | \
    ((val & (1UL << 27)) >> 22) | \
    ((val & (1UL << 16)) >> 10))

  #define HASHKEY8(val) ( \
    ((val & (1UL << 11)) >> 11) | \
    ((val & (1UL << 8)) >> 7) | \
    ((val & (1UL << 14)) >> 12) | \
    ((val & (1UL << 10)) >> 7) | \
    ((val & (1UL << 9)) >> 5) | \
    ((val & (1UL << 27)) >> 22) | \
    ((val & (1UL << 16)) >> 10) | \
    ((val & (1UL << 12)) >> 5))

  // Computes the hash key of a tuple with its dbtuple address.
  static inline uint64_t HashKey(const void *tuple) {
    const uint64_t addrval = reinterpret_cast<uint64_t>(tuple);
    //return ((addrval >> 9) % 1900813) % SetSize;
//    if (BucketsSize == 32) {
//      return HASHKEY5(addrval);
//    } else if (BucketsSize == 64) {
//      return HASHKEY6(addrval);
//    } else if (BucketsSize == 128) {
//      return HASHKEY7(addrval);
//    } else if (BucketsSize == 256) {
//      return HASHKEY8(addrval);
//    } else {
      return (addrval >> 8) & (BucketsSize - 1);
//    }
  }

private:
  size_t n_;
  TupleSetElem tuple_set_[SetSize];
  short hash_first_[BucketsSize];
#ifdef BCC_SHHT_STATS
  size_t n_buckets_occupied_;
#endif
};

template <size_t SetSize>
class VectorTupleSet : public TupleSearchInterface {
public:
  VectorTupleSet() : n_(0) {}

  virtual bool SearchTuple(const void **tuples, int ntuples) const OVERRIDE {
#ifdef BCC_SEARCH_TUPLE_SSE
    ntuples = (ntuples + 1) & (~1);
    for (size_t i = 0; i < n_; ++i) {
      __m128i t1 = _mm_set1_epi64x((uint64_t)(tuple_set_[i]));
      for (int it = 0; it < ntuples; it += 2) {
        __m128i t2 = _mm_load_si128((const __m128i *)&tuples[it]);
        if (_mm_movemask_epi8(_mm_cmpeq_epi64(t1, t2)) != 0)
          return true;
      }
    }
#else
    for (size_t i = 0; i < n_; ++i) {
      for (int it = 0; it < ntuples; ++it) {
        if (tuple_set_[i] == tuples[it])
          return true;
      }
    }
#endif
    return false;
  }

  inline void AddTuple(const void *tuple) {
#ifdef BCC_CHECK
    ALWAYS_ASSERT(n_ < SetSize);
#endif
    tuple_set_[n_++] = tuple;
  }

private:
  size_t n_;
  const void *tuple_set_[SetSize];
};

#endif
