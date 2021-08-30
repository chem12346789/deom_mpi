#include "deom.hpp"
#include "mr.hpp"

void od_corr(DEOMAUX *daux, DEOM *d, const CTRL *c, const char *frho_name, const char *fpolar_name, const char lcr) {
  vector<int> nsave(7);
  FILE *frho = fopen(frho_name, "w");
  FILE *fpolar = fopen(fpolar_name, "w");

  MatrixXcd tmp, tmp1;
  tmp = MatrixXcd::Zero(NSYS, NSYS);
  tmp1 = MatrixXcd::Zero(NSYS, NSYS);
  const double dt = c->dt;
  int nt_i, nt_f, nt;
  double ti;
  if (c->ti < 0) {
    nt_i = ceil(abs(c->ti) / dt);
    nt_f = ceil(abs(c->tf) / dt);
    nt = nt_i + nt_f;
    ti = -nt_i * dt;
  } else {
    nt_i = ceil(abs(c->ti) / dt);
    nt_f = ceil(abs(c->tf) / dt);
    nt = nt_f - nt_i;
    ti = nt_i * dt;
  }

  int i;
  int iado_;
  int re_tree = 1;
  int filter = 1;
  ullint tree_size;

  const int n_threads = omp_get_max_threads();

  printf("===   BRIEFING !!!!   ===\n");
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("tf:%f\n", c->tf);
  printf("n_threads:%d\n", n_threads);

  print(d);

  if (d->dipole.on) {
    print(d->dipole, d);
  } else {
    printf("no dipole, exiting\n");
    exit(122);
  }

  if (d->dipole1.on) {
    print(d->dipole1, d);
  } else {
    printf("no dipole1, exiting\n");
    exit(123);
  }

  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>\n");
  printf("===   BRIEFING DONE   ===\n");

  Trie *tree = new Trie(10 * d->nind);
  tree->try_insert(1, 0);
  const double dt2 = dt / 2.0;
  const double dt3 = dt / 3.0;
  const double dt4 = dt / 4.0;
  const double dt6 = dt / 6.0;

  d->hamt = d->ham1;
  for (int mp = 0; mp < d->nind; mp++) {
    int m1 = d->modLabel(mp);
    d->qmdta[mp].noalias() = d->qmd1a[m1];
    d->qmdtc[mp].noalias() = d->qmd1c[m1];
  }

#pragma omp parallel default(shared)
  {
    // #pragma omp for private(iado_) schedule(dynamic, 64)
    //     for (iado_ = 0; iado_ < d->nddo; iado_++) {
    //       d->g_index_list_pos.row(iado_).setZero();
    //       daux->ddos1[iado_].setZero();
    //     }

#pragma omp single
    {
      nsave[6] = d->nddo;
      tree->clear();
      delete tree;
      Trie *tree = new Trie(min((ullint)2 * (ullint)d->nind * (ullint)d->nddo, d->nmax));
      d->nddo = 0;
      nsave[0] = 0;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[6]; i++)
      filter_p(d, d->ddos, d->keys, daux->ddos1, daux->keys1, tree, i);

#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < nsave[6]; iado_++) {
      d->keys.row(iado_) = daux->keys1.row(iado_);
      d->ddos[iado_].noalias() = daux->ddos1[iado_];
    }

#pragma omp single
    {
      nsave[1] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[1]; i++)
      construct_Mapping_p(d, tree, i);

#pragma omp single
    {
      nsave[2] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[2]; i++)
      rem_oprt(daux->ddos1, d->ddos, d, i, lcr);
#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < nsave[2]; iado_++) {
      d->ddos[iado_].noalias() = daux->ddos1[iado_];
    }
  }

#pragma omp parallel default(shared)
  for (int ii = 0; ii < nt; ii++) {
#pragma omp single
    {
      fprintf(frho, "%16.6e\t", ti + ii * dt);
      int pos = tree->find(1);
      tmp.setZero();

      for (int nmod = 0; nmod < d->nmodmax; nmod++) {
        tmp += d->dipole1.sdip[nmod] * d->ddos[pos];
      }

      // for (int mp = 0; mp < d->nind; mp++) {
      //   int m1 = d->modLabel(mp);

      //   // gen_key(d->zerokey, daux->key_tmp, mp, 1, d);
      //   // pos = tree->find(hash_bad(daux->key_tmp, d));
      //   // if (pos != -1) {
      //   //   int m = d->nmodmaxLabel(mp);
      //   //   if (lcr == 'l') {
      //   //     tmp += (complex<double>)1i * (d->coef_lft[mp](m, ii) * d->dipole.pdipa[mp] * d->ddos[pos]);
      //   //   } else if (lcr == 'r') {
      //   //     tmp += -(complex<double>)1i * (d->coef_rht[mp](m, ii) * d->ddos[pos] * d->dipole.pdipa[mp]);
      //   //   } else if (lcr == 'c') {
      //   //     tmp += (complex<double>)1i * (d->coef_lft[mp](m, ii) * d->dipole.pdipa[mp] * d->ddos[pos] + d->coef_rht[mp](m, ii) * d->ddos[pos] * d->dipole.pdipa[mp]);
      //   //   }
      //   // }

      //   gen_key(d->zerokey, daux->key_tmp, mp, 1, d);
      //   pos = tree->find(hash_bad(daux->key_tmp, d));
      //   if (pos != -1) {
      //     tmp += d->dipole1.bdip(mp) * (complex<double>)1i * (d->dipole1.pdipc[mp] * d->ddos[pos]);
      //   }
      // }

      complex<double> trace =  tmp.trace();

      fprintf(fpolar, "%16.6e", ti + ii * dt);
      fprintf(fpolar, "%20.10e%20.10e", real(trace), imag(trace));
      fprintf(fpolar, "\n");
      pos = tree->find(1);

      for (int i = 0; i < NSYS; i++)
        fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());

      printf("%20.10e%20.10e", real(trace), imag(trace));

      if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
        printf("NH\t");
      } else {
        printf("HE\t");
      }
      // printf("%16.6e\t", d->ddos[pos](1, 0).real());
      // printf("%16.6e\t", d->ddos[pos](1, 0).imag());

      printf("nddo:%d\t", nsave[1]);
      printf("nddo:%d\t", nsave[2]);
      printf("nddo:%d\t", nsave[3]);
      printf("nddo:%d\t", nsave[4]);
      printf("nddo:%d\n", nsave[5]);
    }
#pragma omp single
    {
      nsave[6] = d->nddo;
      ullint tree_size_n = tree->size();
      if (ii < re_tree) {
        tree->clear();
        delete tree;
        tree_size = min((ullint)d->nind * (ullint)d->nind * (ullint)d->nddo, d->nmax);
        tree = new Trie(tree_size);
      } else {
        if (ii % re_tree == 0) {
          tree->clear();
          delete tree;
          tree_size = min((ullint)d->nddo + (ullint)(1.5 * d->nddo), (ullint)2 * d->nmax);
          tree = new Trie(tree_size);
        }
      }

      printf("set tree size: %llu\t real tree size: %llu\t", tree_size, tree_size_n);
      if (ii < filter || ii % filter == 0) {
        d->nddo = 0;
        nsave[0] = 0;
      }
    }

    if (ii < filter || ii % filter == 0) {
#pragma omp for private(iado_) schedule(dynamic, 64)
      for (iado_ = 0; iado_ < nsave[6]; iado_++) {
        d->g_index_list_pos.row(iado_).setZero();
        daux->ddos1[iado_].setZero();
      }

#pragma omp for private(i) schedule(dynamic, 16)
      for (i = 0; i < nsave[6]; i++)
        filter_p(d, d->ddos, d->keys, daux->ddos1, daux->keys1, tree, i);

#pragma omp for private(iado_) schedule(dynamic, 64)
      for (iado_ = 0; iado_ < nsave[6]; iado_++) {
        d->keys.row(iado_) = daux->keys1.row(iado_);
        d->ddos[iado_].noalias() = daux->ddos1[iado_];
      }

#pragma omp single
      nsave[1] = d->nddo;

#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nsave[0]; i < nsave[1]; i++)
        construct_Mapping_p(d, tree, i);

#pragma omp single
      nsave[2] = d->nddo;
    } else {
#pragma omp single
      {
        nsave[1] = nsave[6];
        nsave[2] = nsave[6];
      }
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[2]; i++)
      rem_cal(daux->ddos1, d->ddos, d, i);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++)
      daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[2]; i++)
      daux->ddos3[i].noalias() = daux->ddos1[i] * dt2;

    if (ii < filter || ii % filter == 0) {
#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nsave[0]; i < nsave[3]; i++)
        construct_Mapping_p(daux->ddos3, d, tree, i);

#pragma omp single
      nsave[3] = d->nddo;
    } else {
#pragma omp single
      nsave[3] = nsave[6];
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[3]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos2[i] * dt2;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[2]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = daux->ddos2[i] * dt2;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[2]; i < nsave[3]; i++) {
      daux->ddos1[i].noalias() = daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = daux->ddos2[i] * dt2;
    }

    if (ii < filter || ii % filter == 0) {
#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nsave[0]; i < nsave[3]; i++)
        construct_Mapping_p(daux->ddos3, d, tree, i);

#pragma omp single
      nsave[4] = d->nddo;
    } else {
#pragma omp single
      nsave[4] = nsave[6];
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[4]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos2[i] * dt;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[3]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = daux->ddos2[i] * dt;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[3]; i < nsave[4]; i++) {
      daux->ddos1[i].noalias() = daux->ddos2[i] * 2.0;
      daux->ddos3[i].noalias() = daux->ddos2[i] * dt;
    }

    if (ii < filter || ii % filter == 0) {
#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nsave[0]; i < nsave[4]; i++)
        construct_Mapping_p(daux->ddos3, d, tree, i);

#pragma omp single
      nsave[5] = d->nddo;
    } else {
#pragma omp single
      nsave[5] = nsave[6];
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[5]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i];
      d->ddos[i].noalias() += daux->ddos1[i] * dt6;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[4]; i++) {
      daux->ddos1[i].noalias() += daux->ddos2[i];
      d->ddos[i].noalias() = daux->ddos1[i] * dt6;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[4]; i < nsave[5]; i++) {
      daux->ddos1[i].noalias() = daux->ddos2[i];
      d->ddos[i].noalias() = daux->ddos1[i] * dt6;
    }
  }
  
  fclose(frho);
  fclose(fpolar);
}
