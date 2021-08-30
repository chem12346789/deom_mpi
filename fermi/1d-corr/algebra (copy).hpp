#pragma once
#include "deom.hpp"
#include "index.hpp"

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

llint gen_key(const XXbvet_r &key0, XXbvet_r &key1, const int pos, const int chg, const DEOM *d) {
  ullint hash_value = 1;
  if (!key0(pos) && chg == 1) {
    key1 = key0;
    key1(pos) = true;
    int sum = 0;
    for (int i = 0; i < d->nind; i++) {
      if (key1(i)) {
        sum += 1;
        hash_value += d->comb_list(sum, i) + d->comb_list(sum, d->nind);
      }
    }
  } else if (key0(pos) && chg == -1) {
    key1 = key0;
    key1(pos) = false;
    int sum = 0;
    for (int i = 0; i < d->nind; i++) {
      if (key1(i)) {
        sum += 1;
        hash_value += d->comb_list(sum, i) + d->comb_list(sum, d->nind);
      }
    }
  } else {
    exit(12);
    key1 = d->zerokey;
  }
  return hash_value;
}

int tier(const XXbvet_r &key, const DEOM *d) {
  int tier = 0;
  for (int i = 0; i < d->nind; i++)
    if (key(i))
      tier += 1;
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

llint hash_bad(const XXbvet_r &key, const DEOM *d) {
  ullint hash_value = 1;
  int sum = 0;
  for (int i = 0; i < d->nind; i++) {
    if (key(i)) {
      sum += 1;
      hash_value += d->comb_list(sum, i) + d->comb_list(sum, d->nind);
    }
  }
  return (llint)hash_value;
}