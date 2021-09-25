#pragma once

#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void rem_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado) {
  int pos;
  MatrixNcd ddos;
  ddos = MatrixNcd::Zero(NSYS, NSYS);
  VectorXd key(d->nind);
  pos = d->g_index_list_pos(iado, 0) - 1;
  key = (d->keys.row(iado)).cast<double>();
  complex<double> gams = key.transpose() * d->expn;

  if (pos != -1) {
    ddos.noalias() = -gams * total[pos];
    ddos.noalias() -= (complex<double>)1i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos(iado, mp + 1) - 1;
    if (pos != -1) {
      int m = d->modLabel(mp);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        ddos.noalias() -= (complex<double>)1i * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->qmdt[mp] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->qmdt[mp]);
      }
    }

    pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos.noalias() -= (complex<double>)1i * sqrt(n + 1) * d->coef_abs(mp) * (d->qmdt[mp] * total[pos] - total[pos] * d->qmdt[mp]);
    }
  }

  ddos1[iado].noalias() = ddos;
}

// void syl_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado, const XXimat_r &keys, const XXimat_r &g_index_list_pos) {
//   int pos;
//   double OMG = 0.00001;
//   MatrixNcd ddos;
//   MatrixNcd ddos_old;
//   MatrixNcd ddos_new;
//   ddos = MatrixNcd::Zero(NSYS, NSYS);
//   VectorXd key(d->nind);
//   pos = g_index_list_pos(iado, 0) - 1;
//   ddos_old = total[pos];
//   key = (keys.row(iado)).cast<double>();
//   complex<double> gams = key.transpose() * d->expn;

//   if (pos != -1) {
//     ddos.noalias() = ((complex<double>)1i * OMG - gams) * total[pos];
//   }

//   for (int mp = 0; mp < d->nind; mp++) {
//     double n = key(mp);
//     pos = g_index_list_pos(iado, mp + 1) - 1;
//     if (pos != -1) {
//       int m = d->modLabel(mp);
//       for (int ii = 0; ii < d->nmodmax; ii++) {
//         ddos.noalias() -= (complex<double>)1i * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->qmdt[mp] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->qmdt[mp]);
//       }
//     }

//     pos = g_index_list_pos(iado, mp + 1 + d->nind) - 1;
//     if (pos != -1) {
//       ddos.noalias() -= (complex<double>)1i * sqrt(n + 1) * d->coef_abs(mp) * (d->qmdt[mp] * total[pos] - total[pos] * d->qmdt[mp]);
//     }
//   }

//   MatrixNcd htmp = d->hamt;
//   MatrixNcd hlft = htmp + OMG * MatrixNcd::Identity();
//   MatrixNcd hrht = htmp;
//   ddos = -(complex<double>)1i * ddos;
//   syl_debug(ddos_new, hlft, hrht, ddos);
//   // double diff = max(vectorise(abs(rho_old - ddos.slice(iddo))));
//   // max_diff = max_diff > diff ? max_diff : diff;
//   ddos1[iado].noalias() = ddos_new;
// }

void syl_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado) {
  int pos;
  double OMG = 7;
  MatrixNcd ddos;
  MatrixNcd ddos_old;
  MatrixNcd ddos_new;
  ddos = MatrixNcd::Zero(NSYS, NSYS);
  VectorXd key(d->nind);
  pos = d->g_index_list_pos(iado, 0) - 1;
  ddos_old = total[pos];
  key = (d->keys.row(iado)).cast<double>();
  complex<double> gams = key.transpose() * d->expn;

  if (pos != -1) {
    ddos.noalias() = OMG * total[pos];
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos(iado, mp + 1) - 1;
    if (pos != -1) {
      int m = d->modLabel(mp);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        ddos.noalias() += (complex<double>)1i * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->qmdt[mp] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->qmdt[mp]);
      }
    }

    pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos.noalias() += (complex<double>)1i * sqrt(n + 1) * d->coef_abs(mp) * (d->qmdt[mp] * total[pos] - total[pos] * d->qmdt[mp]);
    }
  }

  MatrixNcd htmp = (complex<double>)1i * d->hamt;
  MatrixNcd hlft = htmp + (OMG + gams) * MatrixNcd::Identity();
  MatrixNcd hrht = htmp;
  syl_complex_debug(ddos_new, hlft, hrht, ddos);
  // double diff = max(vectorise(abs(rho_old - ddos.slice(iddo))));
  // max_diff = max_diff > diff ? max_diff : diff;
  ddos1[iado].noalias() = ddos_new;
}

void rem_oprt(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado, const char lcr) {
  int pos;
  MatrixNcd ddos = MatrixNcd::Zero(NSYS, NSYS);
  XXuivet_r key(d->nind);
  pos = d->g_index_list_pos(iado, 0) - 1;
  key = d->keys.row(iado);

  if (pos != -1) {
    for (int mp = 0; mp < d->nmodmax; mp++) {
      if (lcr == 'l') {
        ddos.noalias() += d->dipole.sdip[mp] * total[pos];
      } else if (lcr == 'r') {
        ddos.noalias() += total[pos] * d->dipole.sdip[mp];
      } else if (lcr == 'c') {
        ddos.noalias() += d->dipole.sdip[mp] * total[pos] - total[pos] * d->dipole.sdip[mp];
      }
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos(iado, mp + 1) - 1;
    if (pos != -1) {
      int m = d->modLabel(mp);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        if (lcr == 'l') {
          ddos.noalias() += d->dipole.bdip[mp] * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->dipole.pdip[m] * total[pos]);
        } else if (lcr == 'r') {
          ddos.noalias() += d->dipole.bdip[mp] * sqrt(n) / d->coef_abs(mp) * (d->coef_rht[mp](m, ii) * total[pos] * d->dipole.pdip[m]);
        } else if (lcr == 'c') {
          ddos.noalias() += d->dipole.bdip[mp] * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->dipole.pdip[m] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->dipole.pdip[m]);
        }
      }
    }

    pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
    if (pos != -1) {
      int m = d->modLabel(mp);
      if (lcr == 'l') {
        ddos.noalias() += d->dipole.bdip[mp] * sqrt(n + 1) * d->coef_abs(mp) * d->dipole.pdip[m] * total[pos];
      } else if (lcr == 'r') {
        ddos.noalias() += d->dipole.bdip[mp] * sqrt(n + 1) * d->coef_abs(mp) * total[pos] * d->dipole.pdip[m];
      } else if (lcr == 'c') {
        ddos.noalias() += d->dipole.bdip[mp] * sqrt(n + 1) * d->coef_abs(mp) * (d->dipole.pdip[m] * total[pos] - total[pos] * d->dipole.pdip[m]);
      }
    }
  }

  ddos1[iado].noalias() = ddos;
}

void construct_Mapping_p(DEOM *d, Trie *tree, int iado) {
  // again and again, I make myself look like a man has the ability to write parallelize code
  // the strategy here is to build a m*pping from iado + mp to right pos
  // first r*duce this question to several small questions that only cal complex<double> number
  // and then label it, the label to the right pos is constructed here
  // note it will not be parallelized due to hight i/o usage(mean it will never be faster than single core), maybe it will be solved in some machine.
  XXuivet_r key(d->nind);
  XXuivet_r key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = d->keys.row(iado);
  hash_val = hash_bad(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos(pos_get, 0) = pos_get + 1;

  if (key.sum() < d->lmax) {
    for (int mp = 0; mp < d->nind; mp++) {
      hash_val = gen_key_p(key, key_swell, mp, 1, d);
      // hash_val = hash_bad(key_swell, d);
      int pos = tree->find(hash_val);
      if (pos == -1) {
        int nddo_tmp;
#pragma omp critical
        {
          nddo_tmp = d->nddo;
          d->nddo++;
        }
        const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
        d->g_index_list_pos(rank_find, mp + 1) = pos_get + 1;
        d->keys.row(rank_find).noalias() = key_swell;
      } else {
        d->g_index_list_pos(pos, mp + 1) = pos_get + 1;
      }
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp) > 0) {
      hash_val = gen_key_n(key, key_swell, mp, -1, d);
      // hash_val = hash_bad(key_swell, d);
      int pos = tree->find(hash_val);
      if (pos == -1) {
        int nddo_tmp;
#pragma omp critical
        {
          nddo_tmp = d->nddo;
          d->nddo++;
        }
        const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
        d->g_index_list_pos(rank_find, mp + d->nind + 1) = pos_get + 1;
        d->keys.row(rank_find).noalias() = key_swell;
      } else {
        d->g_index_list_pos(pos, mp + d->nind + 1) = pos_get + 1;
      }
    }
  }
}

void construct_Mapping_p(nNNmat &ddos, DEOM *d, Trie *tree, int iado) {
  // again and again, I make myself look like a man has the ability to write parallelize code
  // the strategy here is to build a m*pping from iado + mp to right pos
  // first r*duce this question to several small questions that only cal complex<double> number
  // and then label it, the label to the right pos is constructed here
  // note it will not be parallelized due to hight i/o usage(mean it will never be faster than single core), maybe it will be solved in some machine.
  XXuivet_r key(d->nind);
  XXuivet_r key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = d->keys.row(iado);
  hash_val = hash_bad(key, d);
  if (is_valid(ddos, iado, d->ferr) || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos(pos_get, 0) = pos_get + 1;

    if (key.sum() < d->lmax) {
      for (int mp = 0; mp < d->nind; mp++) {
        hash_val = gen_key_p(key, key_swell, mp, 1, d);
        // hash_val = hash_bad(key_swell, d);
        int pos = tree->find(hash_val);
        if (pos == -1) {
          int nddo_tmp;
#pragma omp critical
          {
            nddo_tmp = d->nddo;
            d->nddo++;
          }
          const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
          d->g_index_list_pos(rank_find, mp + 1) = pos_get + 1;
          d->keys.row(rank_find).noalias() = key_swell;
        } else {
          d->g_index_list_pos(pos, mp + 1) = pos_get + 1;
        }
      }
    }

    for (int mp = 0; mp < d->nind; mp++) {
      if (key(mp) > 0) {
        hash_val = gen_key_n(key, key_swell, mp, -1, d);
        // hash_val = hash_bad(key_swell, d);
        int pos = tree->find(hash_val);
        if (pos == -1) {
          int nddo_tmp;
#pragma omp critical
          {
            nddo_tmp = d->nddo;
            d->nddo++;
          }
          const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
          d->g_index_list_pos(rank_find, mp + d->nind + 1) = pos_get + 1;
          d->keys.row(rank_find).noalias() = key_swell;
        } else {
          d->g_index_list_pos(pos, mp + d->nind + 1) = pos_get + 1;
        }
      }
    }
  }
}

void filter_p(DEOM *d, const nNNmat &ddos, const XXuimat_r &keys, nNNmat &ddos1, XXuimat_r &keys1, Trie *tree, int iado) {
  XXuivet_r key(d->nind);
  if (is_valid(ddos[iado], d->ferr)) {
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
