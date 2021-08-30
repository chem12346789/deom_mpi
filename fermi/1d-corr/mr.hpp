// note this part of comment were contained in Principle part of the instruction.(you will get a prefect pdf version)
// Take the boson, for example (Fermion will be written in this form, upto a minus sign at most)
// dotrho^{(n)}_{bf n} = (- i mathcal{L} - sum_{j}n_{j}gamma_{j})rho^{(n)}_{bf n} - i sum_j mathcal{A}_{j} rho^{(n + 1)}_{{bf n}_j^+}- i sum_{j} mathcal{C}_{j} rho^{ (n - 1)}_{{bf n}_j^-},
// The technical details of the super-operator can be determined separately, and there are three things we need to do:
// 1. Store all density matrices
// 2. Look for the density matrix used
// 3. Calculate the super operator and calculate the $dotrho'(n)'/bfn$$
// We first determine a data structure that is mapped by an indicator(e.g.$6_12_20_30_40_5$) to a density matrix label(specifically the location of the physical structure where the density matrix is stored), which we certainly hope will correspond to one by one, or at least one indicator to a matrix label. Because when the indicator is large we need a data structure that is both memory--friendly(red and black) and fast(dictionary tree), choose perfect hash mapping (perfect hash map, PHM), which requires a perfect hash function(Perfect Hash Function, PHF)(strictly defined data structures that are required in parallel conditions, otherwise we only need hash functions that collide as less as possible. Because perfect hash mapping is inherently avoidable racing condition. Let's explain the construction of this function: The first is the construction of the boson :
// The function is selected as $$ \sum_{j = 0} ^ { n - 1 } \binom{j + 1} {sum_{i = 1} ^ { j } n_i + j} $$
// For example, $0 - 11 - 20 - 3 $is $ \binom{1} {0}, \binom{2}{2} \binom{3}{3} = 20$, and the zero indicator is the system indicator $0_10_20_3$, is $\binom{1}{0} \binom{2}{1} \binom{3}{2} = 0$$. It is noted here that we have the risk of having an empty indicator, because if the indicator does not exist it will also be calculated as 0, in the program we calculate the hash of the empty indicator as 0 (and force the program to quit when an empty indicator is encountered). To sum up, 0 is reserved for an empty array, while the hash of the other indicators is added to 1. Fermion is similar in structure.First define $x \otimes n$ as a exclusive or operation, and the function is selected as :
// \sum_{j = 0} ^ { k - 1 } Biggl{binom{sum_{i = 0} ^ j n_j otimes 0} {j} n_j + binom{sum_{i = 0} ^ { j } n_i otimes 0 - 1} {k} n_jBiggr}
// Note that we do not structure Fermi into binary arrays, which are inconvenient to be perfect hash functions at arbitrary truncations, and have obvious and severe hash collision behavior(no matter what the overflow sense).After all the functions have been built successfully, we have the ability to quickly find the upper and lower levels of any indicator so that we can start building parallel architectures.
// To begin by describing the data structure used, we used the following data structures: a hash table that supports atomic operation (provided by Facebook's folly library open source), and array of the storage matrix (provided by eigen library open source) (provided by open source of the c++ standard library).
// Our parallel architecture takes a similar architecture to the mapping protocol, and we begin with the overall process: first we need to build all the information about the upper and lower layers used by any layer metric (e.g. we need $1_12_20_30_40_5$，$0_13_20_30_40_5$，$0_12_21_30_40_5$，$0_12_20_31_40_5$，$0_12_20_30_41_5$and $0_11_20_30_40_5$when calculating). Let's start with the first layer, calculate the hash value of these indicators as a key, and take the corresponding array indecs as a value, use the hash table to find, when the check can not find an array subscript, note that this step is not atomic Operation, of course, is also the part of the entire parallel skeleton that cannot be atomicized, although this step can be made closer to atomic operation by more granular control (each thread uses a different increase cardinality, for example, there are eight threads, then the first thread can only increase by 0). ，8，16。 The second thread can only increase by 1,9,17), but this can result in severe memory consumption so it is not adopted. The resulting array subscript is then stored for backup (index_table). At this point the mapping step is complete, and here is the protocol step.
// Using the index_table obtained from the previous operation, we can calculate the $dot_rho'(n)','bf'$for each layer, because the information on the upper and lower layers is known. The advantage of this is that we avoid the racing condition that we encounter when calculating the sum of the various parts of the $\dot\rho^{(n)}_{\bf n}$."
// Now let us do some translate, in this code g_index_list_pos store all the information about this tier and next or previous tier.
#pragma once

#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"
// rem_cal do reduce calculation of $dot \rho$. Note all the information about this step were contained in ddos(all density matrices) g_index_list_pos(index of next or previous density matrices).
// note the final form of the heom formulas were contained in Principle part of the instruction.
void rem_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado) {
  int pos;
  MatrixNcd ddos;
  ddos = MatrixNcd::Zero(NSYS, NSYS);
  XXbvet_r key(d->nind);
  pos = d->g_index_list_pos(iado, 0) - 1;
  key = d->keys.row(iado);
  complex<double> gams = 0;

  XXbvet_r parity_key(d->nind);
  bool parity = true;
  double parity_p_n;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      gams += d->expn(mp);
      parity ^= true;
    }
    parity_key(mp) = parity;
  }

  if (parity) {
    parity_p_n = 1;
  } else {
    parity_p_n = -1;
  }

  if (pos != -1) {
    ddos.noalias() = -gams * total[pos];
    ddos.noalias() -= (complex<double>)1i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      pos = d->g_index_list_pos(iado, mp + 1) - 1;
      if (pos != -1) {
        int m = d->nmodmaxLabel(mp);
        for (int ii = 0; ii < d->nmodmax; ii++) {
          ddos.noalias() -= (double)(parity_key(mp) ? 1 : -1) * (complex<double>)1i * ((double)parity_p_n * d->coef_lft[mp](m, ii) * d->qmdta[mp] * total[pos] + d->coef_rht[mp](m, ii) * total[pos] * d->qmdta[mp]);
        }
      }
    } else {
      pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
      if (pos != -1) {
        ddos.noalias() -= (double)(parity_key(mp) ? 1 : -1) * (complex<double>)1i * ((double)parity_p_n * d->qmdtc[mp] * total[pos] - total[pos] * d->qmdtc[mp]);
      }
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

// void syl_cal(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado) {
//   int pos;
//   double OMG = 7;
//   MatrixNcd ddos;
//   MatrixNcd ddos_old;
//   MatrixNcd ddos_new;
//   ddos = MatrixNcd::Zero(NSYS, NSYS);
//   VectorXd key(d->nind);
//   pos = d->g_index_list_pos(iado, 0) - 1;
//   ddos_old = total[pos];
//   key = (d->keys.row(iado)).cast<double>();
//   complex<double> gams = key.transpose() * d->expn;

//   if (pos != -1) {
//     ddos.noalias() = OMG * total[pos];
//   }

//   for (int mp = 0; mp < d->nind; mp++) {
//     double n = key(mp);
//     pos = d->g_index_list_pos(iado, mp + 1) - 1;
//     if (pos != -1) {
//       int m = d->modLabel(mp);
//       for (int ii = 0; ii < d->nmodmax; ii++) {
//         ddos.noalias() += (complex<double>)1i * sqrt(n) / d->coef_abs(mp) * (d->coef_lft[mp](m, ii) * d->qmdt[mp] * total[pos] - d->coef_rht[mp](m, ii) * total[pos] * d->qmdt[mp]);
//       }
//     }

//     pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
//     if (pos != -1) {
//       ddos.noalias() += (complex<double>)1i * sqrt(n + 1) * d->coef_abs(mp) * (d->qmdt[mp] * total[pos] - total[pos] * d->qmdt[mp]);
//     }
//   }

//   MatrixNcd htmp = (complex<double>)1i * d->hamt;
//   MatrixNcd hlft = htmp + (OMG + gams) * MatrixNcd::Identity();
//   MatrixNcd hrht = htmp;
//   syl_complex_debug(ddos_new, hlft, hrht, ddos);
//   // double diff = max(vectorise(abs(rho_old - ddos.slice(iddo))));
//   // max_diff = max_diff > diff ? max_diff : diff;
//   ddos1[iado].noalias() = ddos_new;
// }

void rem_oprt(nNNmat &ddos1, const nNNmat &total, const DEOM *d, const int iado, const char lcr) {
  int pos;
  MatrixNcd ddos = MatrixNcd::Zero(NSYS, NSYS);
  XXbvet_r key(d->nind);
  pos = d->g_index_list_pos(iado, 0) - 1;
  key = d->keys.row(iado);

  XXbvet_r parity_key(d->nind);
  bool parity = true;
  double parity_p_n;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      parity ^= true;
    }
    parity_key(mp) = parity;
  }

  if (parity) {
    parity_p_n = 1;
  } else {
    parity_p_n = -1;
  }

  if (pos != -1) {
    for (int mp = 0; mp < d->nmodmax; mp++) {
      if (lcr == 'l') {
        ddos.noalias() += d->dipole.sdip[mp] * total[pos];
      } else if (lcr == 'r') {
        ddos.noalias() += total[pos] * d->dipole.sdip[mp];
      } else if (lcr == 'c') {
        ddos.noalias() += d->dipole.sdip[mp] * total[pos] + total[pos] * d->dipole.sdip[mp];
      }
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      pos = d->g_index_list_pos(iado, mp + 1) - 1;
      if (pos != -1) {
        int m = d->nmodmaxLabel(mp);
        for (int ii = 0; ii < d->nmodmax; ii++) {
          if (lcr == 'l') {
            ddos.noalias() += (double)(parity_key(mp) ? 1 : -1) * (parity_p_n * d->coef_lft[mp](m, ii) * d->dipole.pdipa[mp] * total[pos]);
          } else if (lcr == 'r') {
            ddos.noalias() += -(double)(parity_key(mp) ? 1 : -1) * (d->coef_rht[mp](m, ii) * total[pos] * d->dipole.pdipa[mp]);
          } else if (lcr == 'c') {
            ddos.noalias() += (double)(parity_key(mp) ? 1 : -1) * (parity_p_n * d->coef_lft[mp](m, ii) * d->dipole.pdipa[mp] * total[pos] + d->coef_rht[mp](m, ii) * total[pos] * d->dipole.pdipa[mp]);
          }
        }
      }
    } else {
      pos = d->g_index_list_pos(iado, mp + 1 + d->nind) - 1;
      if (pos != -1) {
        if (lcr == 'l') {
          ddos.noalias() += (double)(parity_key(mp) ? 1 : -1) * (parity_p_n * d->dipole.pdipc[mp] * total[pos]);
        } else if (lcr == 'r') {
          ddos.noalias() += (double)(parity_key(mp) ? 1 : -1) * (total[pos] * d->dipole.pdipc[mp]);
        } else if (lcr == 'c') {
          ddos.noalias() += (double)(parity_key(mp) ? 1 : -1) * (parity_p_n * d->dipole.pdipc[mp] * total[pos] - total[pos] * d->dipole.pdipc[mp]);
        }
      }
    }
  }

  ddos1[iado].noalias() = ddos;
}

void construct_Mapping_p(DEOM *d, Trie *tree, int iado) {
  // construct Mapping and then reduce it, notice that we do different tolerant of cutoff here, because in rk4 algrorithm we need just dt*ddos3/2/1 is valid.
  XXbvet_r key(d->nind);
  XXbvet_r key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = d->keys.row(iado);
  hash_val = hash_bad(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos(pos_get, 0) = pos_get + 1;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      hash_val = gen_key(key, key_swell, mp, -1, d);
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
    } else {
      if (tier(key, d) < d->lmax) {
        hash_val = gen_key(key, key_swell, mp, 1, d);
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
  }
}

void construct_Mapping_p_no_add(DEOM *d, Trie *tree, int iado) {
  // construct Mapping and then reduce it, notice that we do different tolerant of cutoff here, because in rk4 algrorithm we need just dt*ddos3/2/1 is valid.
  XXbvet_r key(d->nind);
  XXbvet_r key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = d->keys.row(iado);
  hash_val = hash_bad(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos(pos_get, 0) = pos_get + 1;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      hash_val = gen_key(key, key_swell, mp, -1, d);
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
  // construct Mapping and then reduce it, notice that we do different tolerant of cutoff here, because in rk4 algrorithm we need just dt*ddos3/2/1 is valid.
  XXbvet_r key(d->nind);
  XXbvet_r key_swell(d->nind);
  llint hash_val;
  // rk4 not rk3 or rk5
  key.noalias() = d->keys.row(iado);
  hash_val = hash_bad(key, d);
  if (is_valid(ddos, iado, d->ferr) || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos(pos_get, 0) = pos_get + 1;
    for (int mp = 0; mp < d->nind; mp++) {
      if (key(mp)) {
        hash_val = gen_key(key, key_swell, mp, -1, d);
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
      } else {
        if (tier(key, d) < d->lmax) {
          hash_val = gen_key(key, key_swell, mp, 1, d);
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
    }
  }
}

void filter_p(DEOM *d, const nNNmat &ddos, const XXbmat_r &keys, nNNmat &ddos1, XXbmat_r &keys1, Trie *tree, int iado) {
  XXbvet_r key(d->nind);
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
