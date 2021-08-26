#include <folly/AtomicHashMap.h>
// #include <map>
// #include <parallel_hashmap/phmap.h>
// #include <sparsehash/dense_hash_map>
// #include <unordered_map>
#include <utility>

using namespace folly;
// using GOOGLE_NAMESPACE::dense_hash_map;
using namespace std;
typedef unsigned long long int ullint;
// using phmap::flat_hash_map;

class Trie {
private:
  // google::dense_hash_map<int64_t, int32_t> tree;
  // unordered_map<int64_t, int32_t> tree;
  // flat_hash_map<int64_t, int32_t> tree;
  AtomicHashMap<int64_t, int32_t> tree;
  // Try to insert a key-value pair; return true if succeed, false if exist;
public:
  explicit Trie(size_t nmax) : tree((ullint)1 * (ullint)nmax) {}
  // explicit Trie(size_t nmax) {}
  // explicit Trie(ullint nmax) { tree.set_empty_key((int64_t)(nmax + 1)); }

  void insert(const int64_t key_hash, const int32_t rank) {
    tree.insert({key_hash, rank});
  }

  pair<int32_t, bool> try_insert(const int64_t key_hash, const int32_t rank) {
    const auto [rank_find, success] = tree.insert({key_hash, rank});
    return make_pair(rank_find->second, success);
  }

  int32_t find(const int64_t key_hash) const {
    auto search = tree.find(key_hash);
    if (search != tree.end())
      return search->second;
    else
      return -1;
  }

  void remove(const int64_t key_hash) {
    auto search = tree.find(key_hash);
    if (search != tree.end()) {
      bool success = tree.erase(key_hash);
      if (!success) {
        printf("ERROR IN MAP!\n");
        exit(1);
      }
    }
  }

  void clear() {
    tree.clear();
  }

  int size() {
    return tree.size();
  }
};