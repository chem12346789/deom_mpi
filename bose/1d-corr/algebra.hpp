#pragma once
#include "deom.hpp"
#include "index.hpp"

// const complex mat, double > bool
// weather the key is valid(all item are greater than ferr(ie. 1e-10))
bool is_valid(const MatrixNcd &d_ddo, const double ferr = 1e-10) {
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      if (abs(d_ddo(i, j)) > ferr)
        return true;
  return false;
}

bool is_valid(const nNNmat &d_ddos, const int iado, const double ferr = 1e-10) {
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      if (abs(d_ddos[iado](i, j)) > ferr)
        return true;
  return false;
}

// const uint array, uint array, const int, const int, const DEOM > hash and const uint array, uint array(gened key), ...
// Haskell done now turn to c++.
// gen the plus anf minus indicator.
// we choose Perfect Hash Function as
// $$
// \sum_{j = 1} ^ { n }\binom{j + 1} {\sum_{i = 1} ^ { j } n_i + j}
// $$
// EXP. $0_11_20_3$ as $\binom{1} {0} +\binom{2} {2} +\binom{3} {3} = 2 $，zero key also system key is $0_10_20_3 = $ as $\binom{1} {0} +\binom{2} {1} +\binom{3} {2} = 0 $。It is noted here that we have the risk of having an empty indicator, because if the indicator does not exist it will also be calculated as 0, in the program we calculate the hash of the indicator as 0 (and force the end of the program when an empty indicator is encountered). To sum up, 0 is reserved for an empty array, while the hash of the other indicators is added to one.
llint gen_key(const XXuivet_r &key0, XXuivet_r &key1, const int pos, const int chg, const DEOM *d) {
  ullint hash_value = 1;
  if (chg > 0) {
    key1 = key0;
    key1(pos) += chg;
    int sum = 0;
    for (int i = 0; i < d->nind; i++) {
      sum += key1(i);
      hash_value += d->comb_list(sum + i, i + 1);
    }
  } else if (key0(pos) + chg >= 0) {
    key1 = key0;
    key1(pos) += chg;
    int sum = 0;
    for (int i = 0; i < d->nind; i++) {
      sum += key1(i);
      hash_value += d->comb_list(sum + i, i + 1);
    }
  } else {
    exit(0);
    // when we meet empty key, we stop, and look back in anger.
  }
  return hash_value;
}

int tier(const XXuivet_r &key, const DEOM *d) {
  int tier = 0;
  for (int i = 0; i < d->nind; i++)
    tier += key(i);
  return tier;
}

int comb(const int x, const int y) {
  if (x >= y) {
    if (x <= 2 * y) {
      ullint res = 1;
      for (int i = y + 1; i <= x; i++)
        res = res * i;
      for (int i = 1; i <= x - y; i++)
        res = res / i;
      return (int)res;
    } else {
      ullint res = 1;
      for (int i = x - y + 1; i <= x; i++)
        res = res * i;
      for (int i = 1; i <= y; i++)
        res = res / i;
      return (int)res;
    }
  } else {
    return 0;
  }
}

llint hash_bad(const XXuivet_r &key, const DEOM *d) {
  ullint hash_value = 1;
  // why 1?
  // because i do not like 0, it was special hash value for empty key.
  int sum = 0;
  for (int i = 0; i < d->nind; i++) {
    sum += key(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return (llint)hash_value;
}