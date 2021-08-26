#include "deom.hpp"
#include "index.cpp"

bool is_valid(const MatrixXcd &d_ddo, const double ferr = 1e-10) {
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

void gen_key(const VectorXi &key0, VectorXi &key1, const int pos, const int chg, const DEOM *d) {
  int npos, len0 = d->nind;

  for (int i = 0; i < d->nind; i++) {
    if (key0(i) == -1) {
      len0 = i;
      break;
    }
  }

  if (len0 >= pos + 1 && chg > 0) {
    for (int i = 0; i < len0; i++)
      key1(i) = key0(i);
    key1(pos) += chg;
    for (int i = len0; i < d->nind; i++)
      key1(i) = -1;
  } else if (len0 < pos + 1 && chg > 0) {
    for (int i = 0; i < len0; i++)
      key1(i) = key0(i);
    for (int i = len0; i < pos + 1; i++)
      key1(i) = 0;
    for (int i = pos + 1; i < d->nind; i++)
      key1(i) = -1;
    key1(pos) += chg;
  } else if (len0 >= pos + 1 && chg < 0 && key0(pos) + chg >= 0) {
    for (int i = 0; i < len0; i++)
      key1(i) = key0(i);
    key1(pos) += chg;
    npos = len0;
    while (key1(npos - 1) == 0 and npos != 1)
      npos -= 1;
    for (int i = npos; i < d->nind; i++)
      key1(i) = -1;
  } else {
    for (int i = 0; i < d->nind; i++)
      key1(i) = -1;
  }
}

int tier(const VectorXi &key, const DEOM *d) {
  int tier = 0;
  for (int i = 0; i < d->nind; i++)
    if (key(i) != -1)
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

llint hash_bad(const VectorXi &key, const DEOM *d) {
  ullint hash_value = 0;
  if (key(0) != -1) {
    hash_value = 1;
  }
  int sum = 0;
  for (int i = 0; i < d->nind; i++) {
    if (key(i) != -1) {
      sum += key(i);
    }
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return (llint)hash_value;
}

