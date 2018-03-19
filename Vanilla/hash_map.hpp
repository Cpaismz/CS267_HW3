#pragma once

#include <vector>
#include <mutex>
#include <upcxx/upcxx.hpp>
#include "kmer_t.hpp"

struct HashMap {
    std::vector<upcxx::global_ptr<kmer_pair>> data;
    std::vector<upcxx::global_ptr<int>> used;

  size_t my_size;
  std::mutex lock;

  size_t size() const noexcept;

  HashMap(size_t size, std::vector<upcxx::global_ptr<kmer_pair>> _data, std::vector<upcxx::global_ptr<int>> _used);

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  void insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

  // Helper functions

  // Write and read to a logical data slot in the table.
  void write_slot(uint64_t slot, const kmer_pair &kmer, int node_num);
  kmer_pair read_slot(uint64_t slot, int node_num);

  // Request a slot or check if it's already used.
  bool request_slot(uint64_t slot, int node_num);
  bool slot_used(uint64_t slot, int node_num);
};

HashMap::HashMap(size_t size, std::vector<upcxx::global_ptr<kmer_pair>> _data, std::vector<upcxx::global_ptr<int>> _used) {
  my_size = size;
  data = _data;
  used = _used;
}


// only insert "locally"
void HashMap::insert(const kmer_pair &kmer) {
  lock.lock();
  uint64_t hash = kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  int node_num = (hash % size()) / (size() / upcxx::rank_n());
  do {
    uint64_t slot = (hash + probe++) % size();
    success = request_slot(slot, node_num);
    if (success) {
      write_slot(slot, kmer, node_num);
    }
  } while (!success && probe < size());
  lock.unlock();
}

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
  uint64_t hash = key_kmer.hash();
  uint64_t probe = 0;
  int node_num = (hash % size()) / (size() / upcxx::rank_n());
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    if (slot_used(slot, node_num)) {
      val_kmer = read_slot(slot, node_num);
      if (val_kmer.kmer == key_kmer) {
        success = true;
      }
    }
  } while (!success && probe < size());
  return success;
}

bool HashMap::slot_used(uint64_t slot, int node_num) {
    return rget(used[node_num] + slot).wait() != 0;
  // return used[slot] != 0;
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer, int node_num) {
    rput(kmer, data[node_num] + slot).wait();
}


kmer_pair HashMap::read_slot(uint64_t slot, int node_num) {
  return rget(data[node_num] + slot).wait();
  //return data[slot];
}

bool HashMap::request_slot(uint64_t slot, int node_num) {
  if (rget(used[node_num] + slot).wait() == 0) {
  // if (used[slot] != 0) {
    return false;
  } else {
    rput(1, used[node_num] + slot).wait();
    // used[slot] = 1;
    return true;
  }
}

size_t HashMap::size() const noexcept {
  return my_size;
}
