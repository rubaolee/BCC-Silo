#ifndef __BCC_MM_H__
#define __BCC_MM_H__

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <malloc.h>

#include <iostream>

#include "list.h"
#include "macros.h"
#include "sh_macros.h"

// Use BCC block management; NOT recommended.
//#define BCC_BLOCK_MM

#ifdef BCC_BLOCK_MM
// Max number of mm blocks. Must be power of 2.
#define BLOCKNUM    (64 * 1024)
#define BLOCKMASK   (BLOCKNUM - 1)
#endif

#ifdef BCC_SHMM_STATS
struct shmm_stats {
  void init() {
#ifdef BCC_BLOCK_MM
    count_blk_alloc = 0;
    count_blk_free = 0;
    max_allocated_blks = 0;
    curr_allocated_blks = 0;
#endif
    time_alloc = 0.0;
    time_release = 0.0;
    total_allocated_mem = 0;
    nallocs = 0;
    max_allocated_mem = 0;
    curr_allocated_mem = 0;
  }

  void print() {
    std::cerr << "SHMM STATS: "
#ifdef BCC_BLOCK_MM
              << "count_blk_alloc(" << count_blk_alloc << ") "
              << "count_blk_free(" << count_blk_free << ") "
              << "max_allocated_blks(" << max_allocated_blks << ") "
#endif
              << "time_alloc(" << time_alloc << ") "
              << "time_release(" << time_release << ") "
              << "avg_allocated_mem("
              << (nallocs > 0 ? (total_allocated_mem / nallocs) : 0) << ") "
              << "max_allocated_mem(" << max_allocated_mem << ") "
              << std::endl;
  }

#ifdef BCC_BLOCK_MM
  void inc_block_alloc() {
    ++count_blk_alloc;
    if (++curr_allocated_blks > max_allocated_blks)
      max_allocated_blks = curr_allocated_blks;
  }

  void inc_block_free() {
    ++count_blk_free;
    --curr_allocated_blks;
  }
#endif

  void add_time_alloc(double delta) {
    time_alloc += delta;
  }

  void add_time_release(double delta) {
    time_release += delta;
  }

  void add_size_alloc(size_t size) {
    curr_allocated_mem += size;
    total_allocated_mem += curr_allocated_mem;
    ++nallocs;
    if (curr_allocated_mem > max_allocated_mem)
      max_allocated_mem = curr_allocated_mem;
  }

  void add_size_release(size_t size) {
    curr_allocated_mem -= size;
  }

#ifdef BCC_BLOCK_MM
  // How many times block_alloc/block_free have been called.
  int64_t count_blk_alloc;
  int64_t count_blk_free;
  // Maximum number of allocated blocks at any time.
  int64_t max_allocated_blks;
  // Current number of allocated blocks; used only for computing
  // max_allocated_blks.
  // *** Do NOT print this field ***
  int64_t curr_allocated_blks;
#endif
  // Total time (us) spent in shmm allocation/release.
  double time_alloc;
  double time_release;
  // For computing average size of allocated memory.
  size_t total_allocated_mem;
  uint64_t nallocs;
  // Maximum size of allocated memory at any time.
  size_t max_allocated_mem;
  // Current size of allocated memory; used only for computing
  // max_allocated_mem.
  // *** Do NOT print this field ***
  size_t curr_allocated_mem;
};
#endif

// NOTE: Do NOT create an object of this class on the stack ifdef BCC_BLOCK_MM!!
class shmm {
protected:
  struct block {
#ifdef BCC_BLOCK_MM
    int index;
#endif
    char *start;
    size_t size;
    struct list_head list;  // Link to list_free or list_alloc
  };

  inline struct block *block_alloc() {
#ifdef BCC_BLOCK_MM
    for (int i = blk_alloc_search_startidx, k = 0; k < BLOCKNUM;
        i = (i + 1) & BLOCKMASK, ++k) {
      if (block_status[i] == 0) {
        block_status[i] = 1;
        block_buf[i].index = i;
        blk_alloc_search_startidx = (i + 1) & BLOCKMASK;
        stats.inc_block_alloc();
        return &block_buf[i];
      }
    }
    std::cout << "block_alloc fail: "
              << stats.count_blk_alloc << " " << stats.count_blk_release
              << std::endl;
    ALWAYS_ASSERT(false);
    return nullptr;
#else
    return (struct block *)malloc(sizeof(struct block));
#endif
  }

  inline void block_free(struct block * blk) {
#ifdef BCC_BLOCK_MM
    block_status[blk->index] = 0;
    stats.inc_block_free();
#else
    free(blk);
#endif
  }

public:
  void init(int worker_idx /* unused */, char *ad, size_t size) {
#ifdef BCC_CHECK
    // Currently, only supports three mempools.
    ALWAYS_ASSERT(SHMM_MEMPOOLS == 3);
    ALWAYS_ASSERT(SHMM_MEMPOOL_SIZE_0 + SHMM_MEMPOOL_SIZE_1 + SHMM_MEMPOOL_SIZE_2 <= size);
#endif
#ifdef BCC_BLOCK_MM
    blk_alloc_search_startidx = 0;
    blk_alloc_count = 0;
    blk_release_count = 0;
#endif
#ifdef BCC_SHMM_STATS
    stats.init();
#endif

    for (int ipool = 0; ipool < SHMM_MEMPOOLS; ++ipool) {
      INIT_LIST_HEAD(&mempool[ipool].list_free);
      INIT_LIST_HEAD(&mempool[ipool].list_alloc);
    }
#ifdef BCC_BLOCK_MM
    memset(block_status, 0, sizeof(block_status));
#endif
    struct block *newblk = block_alloc();
    newblk->start = ad;
    newblk->size = SHMM_MEMPOOL_SIZE_0;
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    last_alloc_end = ad;
#endif
    list_add(&newblk->list, &mempool[0].list_free);
    // Warm up cache for mempool 0
//    int64_t dummy = 0;
//    for (char *p = ad; p < ad + SHMM_MEMPOOL_SIZE_0; p += CACHELINE_SIZE)
//      dummy += *p;
//    *(int64_t *)ad = dummy;
    newblk = block_alloc();
    newblk->start = ad + SHMM_MEMPOOL_SIZE_0;
    newblk->size = SHMM_MEMPOOL_SIZE_1;
    list_add(&newblk->list, &mempool[1].list_free);
    newblk = block_alloc();
    newblk->start = ad + SHMM_MEMPOOL_SIZE_0 + SHMM_MEMPOOL_SIZE_1;
    newblk->size = SHMM_MEMPOOL_SIZE_2;
    list_add(&newblk->list, &mempool[2].list_free);
    memset(&inflight, 0, sizeof(inflight));
  }

#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
  char *volatile last_alloc_end;
#endif
  /*
   * First try to allocate from mempool[0], and, if failed, then try to
   * allocate from mempool[1].
   */
  std::pair<char *, int> allocate(size_t size) {
    std::pair<char *, int> ret(nullptr, 0);
#ifdef BCC_DEBUG_DISABLE_SHMM_SHTXN_OPS
    struct block *b = list_entry(mempool[0].list_free.next, struct block, list);
    if (last_alloc_end + size > b->start + b->size)
      last_alloc_end = b->start;
    ret.first = last_alloc_end;
    last_alloc_end += size;
    ret.second = 0;
    return ret;
#endif
    for (int ipool = 0; ipool < SHMM_MEMPOOLS; ++ipool) {
      ret.first = mempool_allocate(ipool, size);
      ret.second = ipool;
      if (ret.first)
        break;
    }
#ifdef BCC_CHECK
//    if (!ret.first) {
//      // This should not happen; print mempools for debugging.
//      for (int ipool = 0; ipool < SHMM_MEMPOOLS; ++ipool) {
//        print_mempool(ipool);
//      }
//    }
    ALWAYS_ASSERT(ret.first);
#endif
#ifdef BCC_SHMM_STATS
    stats.add_size_alloc(size);
#endif
    return ret;
  }

  char *mempool_allocate(int ipool, size_t size) {
    struct list_head &list_free = mempool[ipool].list_free;
    struct list_head &list_alloc = mempool[ipool].list_alloc;
    char *ret = nullptr;
#ifdef BCC_SHMM_STATS
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
#endif

    // Search free list from head to tail.
    struct list_head *pos = list_free.next;
    while (pos != &list_free) {
      struct block *b = list_entry(pos, struct block, list);
      if (likely(b->size > size)) {
        set_inflight(b->start, size);
        ret = b->start;
        b->start += size;
        b->size -= size;
        add_inflight_into_list(&list_alloc);
        break;
      } else if (unlikely(b->size == size)) {
        list_del(&b->list);
        set_inflight(b);
        ret = b->start;
        add_inflight_into_list(&list_alloc);
        break;
      } else {
        // When a free block becomes too small, move it to list_alloc.
#ifdef BCC_SHMM_STATS
        stats.add_size_alloc(b->size);
#endif
        pos = pos->next;
        list_del(&b->list);
        set_inflight(b);
        add_inflight_into_list(&list_alloc);
      }
    }

#ifdef BCC_SHMM_STATS
    clock_gettime(CLOCK_REALTIME, &tend);
    stats.add_time_alloc(TDIFF(tstart, tend));
#endif
    return ret;
  }

  /*
   * Release allocated memory until @end_addr.
   */
  void release(int ipool, const char *end_addr) {
#ifdef BCC_CHECK
    ALWAYS_ASSERT(ipool >= 0 && ipool < SHMM_MEMPOOLS);
#endif
    struct list_head &list_free = mempool[ipool].list_free;
    struct list_head &list_alloc = mempool[ipool].list_alloc;
#ifdef BCC_SHMM_STATS
    size_t size_released = 0;
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
#endif

    struct list_head *pos = list_alloc.next;
    while (pos != &list_alloc) {
      struct block *b = list_entry(pos, struct block, list);
      if (end_addr > b->start && end_addr < (b->start + b->size)) {
        // Part of this allocated block is to be released.
        size_t size_delta = (size_t)end_addr - (size_t)(b->start);
#ifdef BCC_SHMM_STATS
        size_released += size_delta;
#endif
        set_inflight(b->start, size_delta);
        b->start += size_delta;
        b->size -= size_delta;
        add_inflight_into_list(&list_free);
        break;
      } else {
        // The whole allocated block is to be released.
//#ifdef BCC_CHECK
//        ALWAYS_ASSERT(end_addr == b->start + b->size || end_addr <= b->start);
//#endif
        pos = pos->next;
        list_del(&b->list);
#ifdef BCC_SHMM_STATS
        size_released += b->size;
#endif
        const bool do_break = (end_addr == b->start + b->size);
        set_inflight(b);
        add_inflight_into_list(&list_free);
        if (do_break)
          break;
      }
    }

#ifdef BCC_SHMM_STATS
    clock_gettime(CLOCK_REALTIME, &tend);
    stats.add_time_release(TDIFF(tstart, tend));
    stats.add_size_release(size_released);
#endif
  }

  void print_stats() {
#ifdef BCC_SHMM_STATS
    stats.print();
#endif
  }

protected:
  bool merge_into_list(
      struct list_head *list_merge, const char *start_addr, size_t size) {
    if (list_empty(list_merge))
      return false;
    struct block *tail_blk = list_entry(list_merge->prev, struct block, list);
    if (tail_blk->start + tail_blk->size == start_addr) {
      tail_blk->size += size;
      return true;
    }
    return false;
  }

  bool merge_inflight_into_list(struct list_head *list_merge) {
    ALWAYS_ASSERT(inflight.start_addr);
    if (merge_into_list(list_merge, inflight.start_addr, inflight.size)) {
      if (inflight.inflight_block) {
        block_free(inflight.inflight_block);
      }
      return true;
    }
    return false;
  }

  void insert_inflight_into_list(struct list_head *list_insert) {
    if (inflight.inflight_block) {
      // Reuse free_block: directly insert it into list_alloc.
      list_add_tail(&inflight.inflight_block->list, list_insert);
    } else {
      // A new block structure needs to be allocated.
      struct block *newblk = block_alloc();
      newblk->start = inflight.start_addr;
      newblk->size = inflight.size;
      list_add_tail(&newblk->list, list_insert);
    }
  }

  /*
   * Add inflight information into @list_add (which could be either
   * list_free or list_alloc). Merge with an existing block on
   * @list_add if possible.
   */
  void add_inflight_into_list(struct list_head *list_add) {
    if (!merge_inflight_into_list(list_add)) {
      insert_inflight_into_list(list_add);
    }
  }

  void set_inflight(struct block *inflight_block) {
#ifdef BCC_CHECK
    ALWAYS_ASSERT(inflight_block);
#endif
    inflight.start_addr = inflight_block->start;
    inflight.size = inflight_block->size;
    inflight.inflight_block = inflight_block;
  }

  void set_inflight(char *start_addr, size_t size) {
    inflight.start_addr = start_addr;
    inflight.size = size;
    inflight.inflight_block = nullptr;
  }

  // For debugging purpose only
  void print_mempool(int ipool) {
    std::cerr << "Printing mempool " << ipool << ":" << std::endl;
    std::cerr << "Free list (head->tail):" << std::endl;
    print_list(&mempool[ipool].list_free);
    std::cerr << "Alloc list (head->tail):" << std::endl;
    print_list(&mempool[ipool].list_alloc);
  }

  // For debugging purpose only
  void print_list(struct list_head *list_print) {
    struct list_head *pos;
    list_for_each(pos, list_print) {
      struct block *b = list_entry(pos, struct block, list);
      std::cerr << "(" << reinterpret_cast<void *>(b->start) << ","
                << reinterpret_cast<void *>(b->start + b->size) << ","
                << b->size << ") ";
    }
    std::cerr << std::endl;
  }

private:
#ifdef BCC_BLOCK_MM
  // Search for unused block starts here
  int blk_alloc_search_startidx;
  struct block block_buf[BLOCKNUM];
  char block_status[BLOCKNUM];
#endif

  // The locality of data managed by mempools decreases as the index
  // of mempool increases. For example, mempool[0]'s locality is stronger
  // than mempool[1]'s, which is stronger than mempool[2]'s, ....
  struct {
    // List of free blocks.
    // The list head always points to the block where search for new allocation
    // should begin. In current implementation, there should be at most two free
    // blocks at any time.
    struct list_head list_free;
    // List of allocated blocks.
    // Allocated blocks are sorted in the increasing order of their epochs.
    // List head points to the block with the smallest epoch; list tail points
    // to the block with the largest. Two adjacent blocks may have the same
    // epoch, in the case that the block start address wraps around the end of
    // managed memory area. Newly allocated blocks should always be inserted
    // or merged to the list tail.
    struct list_head list_alloc;
  } mempool[SHMM_MEMPOOLS];
  // Stores information about a block of memory that has been allocated but
  // not been inserted to list_alloc yet.
  struct {
    char *start_addr;
    size_t size;
    // If inflight_block is not nullptr, (start_addr, size) is directly taken
    // from the block and inflight owns inflight_block in this case.
    struct block *inflight_block;
  } inflight;

#ifdef BCC_SHMM_STATS
  shmm_stats stats;
#endif
};

#endif
