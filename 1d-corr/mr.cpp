#include "algebra.cpp"
#include "deom.hpp"

void rem_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado, const XXimat_r &keys, const XXimat_r &g_index_list_pos) {
  int pos;
  MatrixXcd ddos;
  ddos = MatrixXcd::Zero(NSYS, NSYS);

  VectorXi key(d->nind);
  complex<double> gams = 0;

  pos = g_index_list_pos(iado, 0) - 1;

  key.noalias() = keys.row(iado);

  for (int k = 0; k < d->nind; k++) {
    if (key[k] == -1)
      break;
    gams += (double)key(k) * d->expn(k);
  }

  if (pos != -1) {
    ddos.noalias() = -gams * total[pos];
    ddos.noalias() -= (complex<double>)1i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp) == -1 ? 0 : key(mp);
    pos = g_index_list_pos(iado, mp + 1) - 1;
    if (pos != -1) {
      int m = d->modLabel(mp);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        ddos.noalias() -= (complex<double>)1i * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->qmdt[mp] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->qmdt[mp]);
      }
    }

    pos = g_index_list_pos(iado, mp + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos.noalias() -= (complex<double>)1i * sqrt(n + 1) * d->coef_abs(mp) * (d->qmdt[mp] * total[pos] - total[pos] * d->qmdt[mp]);
    }
  }

  ddos1[iado].noalias() = ddos;
}

void construct_Mapping_p(DEOM *d, XXimat_r &keys, XXimat_r &g_index_list_pos, Trie *tree, int iado) {
  // again and again, I make myself look like a man has the ability to write parallelize code
  // the strategy here is to build a m*pping from iado + mp to right pos
  // first r*duce this question to several small questions that only cal complex<double> number
  // and then label it, the label to the right pos is constructed here
  // note it will not be parallelized due to hight i/o usage(mean it will never be faster than single core), maybe it will be solved in some machine.
  VectorXi key(d->nind);
  VectorXi key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = keys.row(iado);
  hash_val = hash_bad(key, d);
  int pos_get = tree->find(hash_val);
  g_index_list_pos(pos_get, 0) = pos_get + 1;

  if (tier(key, d) < d->lmax) {
    for (int mp = 0; mp < d->nind; mp++) {
      gen_key(key, key_swell, mp, 1, d);
      hash_val = hash_bad(key_swell, d);
      int pos = tree->find(hash_val);
      if (pos == -1) {
        int nddo_tmp;
#pragma omp critical
        {
          nddo_tmp = d->nddo;
          d->nddo++;
        }
        const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
        g_index_list_pos(rank_find, mp + 1) = pos_get + 1;
        keys.row(rank_find).noalias() = key_swell;
      } else {
        g_index_list_pos(pos, mp + 1) = pos_get + 1;
      }
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp) > 0) {
      gen_key(key, key_swell, mp, -1, d);
      hash_val = hash_bad(key_swell, d);
      int pos = tree->find(hash_val);
      if (pos == -1) {
        int nddo_tmp;
#pragma omp critical
        {
          nddo_tmp = d->nddo;
          d->nddo++;
        }
        const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
        g_index_list_pos(rank_find, mp + d->nind + 1) = pos_get + 1;
        keys.row(rank_find).noalias() = key_swell;
      } else {
        g_index_list_pos(pos, mp + d->nind + 1) = pos_get + 1;
      }
    }
  }
}

void construct_Mapping_p(nNNmat &ddos, DEOM *d, XXimat_r &keys, XXimat_r &g_index_list_pos, Trie *tree, int iado) {
  // again and again, I make myself look like a man has the ability to write parallelize code
  // the strategy here is to build a m*pping from iado + mp to right pos
  // first r*duce this question to several small questions that only cal complex<double> number
  // and then label it, the label to the right pos is constructed here
  // note it will not be parallelized due to hight i/o usage(mean it will never be faster than single core), maybe it will be solved in some machine.
  VectorXi key(d->nind);
  VectorXi key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = keys.row(iado);
  hash_val = hash_bad(key, d);
  if (is_valid(ddos, iado, d->ferr) || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    g_index_list_pos(pos_get, 0) = pos_get + 1;

    if (tier(key, d) < d->lmax) {
      for (int mp = 0; mp < d->nind; mp++) {
        gen_key(key, key_swell, mp, 1, d);
        hash_val = hash_bad(key_swell, d);
        int pos = tree->find(hash_val);
        if (pos == -1) {
          int nddo_tmp;
#pragma omp critical
          {
            nddo_tmp = d->nddo;
            d->nddo++;
          }
          const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
          g_index_list_pos(rank_find, mp + 1) = pos_get + 1;
          keys.row(rank_find).noalias() = key_swell;
        } else {
          g_index_list_pos(pos, mp + 1) = pos_get + 1;
        }
      }
    }

    for (int mp = 0; mp < d->nind; mp++) {
      if (key(mp) > 0) {
        gen_key(key, key_swell, mp, -1, d);
        hash_val = hash_bad(key_swell, d);
        int pos = tree->find(hash_val);
        if (pos == -1) {
          int nddo_tmp;
#pragma omp critical
          {
            nddo_tmp = d->nddo;
            d->nddo++;
          }
          const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
          g_index_list_pos(rank_find, mp + d->nind + 1) = pos_get + 1;
          keys.row(rank_find).noalias() = key_swell;
        } else {
          g_index_list_pos(pos, mp + d->nind + 1) = pos_get + 1;
        }
      }
    }
  }
}

void filter_p(DEOM *d, const nNNmat &ddos, const XXimat_r &keys, nNNmat &ddos1, XXimat_r &keys1, Trie *tree, int iado) {
  VectorXi key(d->nind);
  if (is_valid(ddos, iado, d->ferr)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    key.noalias() = keys.row(iado);
    llint hash_val = hash_bad(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp].noalias() = ddos[iado];
    keys1.row(nddo_tmp).noalias() = key;
  }
}
