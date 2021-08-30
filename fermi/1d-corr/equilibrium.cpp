#include "deom.hpp"
#include "mr.hpp"

void equilibrium(DEOMAUX *daux, DEOM *d, CTRL *c, const char *frho_name) {
  vector<int> nsave(7);
  FILE *frho = fopen(frho_name, "w");
  FILE *fcur = fopen("curr.dat", "w");

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

  int re_tree = 10;
  int filter = 10;

  const int n_threads = omp_get_max_threads();
  ullint tree_size;

  printf("===   BRIEFING !!!!   ===\n");
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("tf:%f\n", c->tf);
  print(d);

  printf("task: equilibrium, cal \\rho^{eq}\n");
  printf("===   BRIEFING DONE   ===\n");

  Trie *tree = new Trie(d->nmax);
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

  // int pos = tree->find(1);
  // d->ddos[pos](0, 0) = 0.162102;
  // d->ddos[pos](1, 1) = 0.337898;
  // d->ddos[pos](2, 2) = 0.337898;
  // d->ddos[pos](3, 3) = 0.162102;

#pragma omp parallel default(shared)
  for (int ii = 0; ii < nt; ii++) {
#pragma omp single
    {
      fprintf(frho, "%16.6e", ti + ii * dt);
      int pos = tree->find(1);
      for (int i = 0; i < NSYS; i++)
        fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());
      fprintf(frho, "\n");
      printf("%16.6e\t", ti + ii * dt);
      for (int i = 0; i < NSYS; i++)
        printf("%16.6e\t", d->ddos[pos](i, i).real());

      fprintf(fcur, "%16.6e\t", ti + ii * dt);
      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        gen_key(d->zerokey, daux->key_tmp, mp, 1, d);
        pos = tree->find(hash_bad(daux->key_tmp, d));
        if (pos != -1) {
          fprintf(fcur, "%16.6e\t%16.6e\t", ((complex<double>)1i * (d->qmdtc[mp] * d->ddos[pos]).trace()).real(), ((complex<double>)1i * (d->qmdtc[mp] * d->ddos[pos]).trace()).imag());
        } else {
          fprintf(fcur, "%16.6e\t%16.6e\t", 0.0, 0.0);
        }
      }
      fprintf(fcur, "\n");

      if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
        printf("NH\t");
      } else {
        printf("HE\t");
      }
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
          tree_size = min((ullint)d->nddo + (ullint)(0.5 * d->nddo), (ullint)2 * d->nmax);
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
      {
        nsave[1] = d->nddo;
        printf("nddo:%12d", nsave[1]);
      }

#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nsave[0]; i < nsave[1]; i++)
        construct_Mapping_p(d, tree, i);

#pragma omp single
      {
        nsave[2] = d->nddo;
        printf("nddo:%12d", nsave[2]);
      }
    } else {
#pragma omp single
      {
        nsave[1] = nsave[6];
        nsave[2] = nsave[6];
        printf("nddo:%12d", nsave[1]);
        printf("nddo:%12d", nsave[2]);
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
      for (i = nsave[0]; i < nsave[2]; i++)
        construct_Mapping_p(daux->ddos3, d, tree, i);

#pragma omp single
      {
        nsave[3] = d->nddo;
        printf("nddo:%12d", nsave[3]);
      }

    } else {
#pragma omp single
      {
        nsave[3] = nsave[6];
        printf("nddo:%12d", nsave[3]);
      }
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
      {
        nsave[4] = d->nddo;
        printf("nddo:%12d", nsave[4]);
      }
    } else {
#pragma omp single
      {
        nsave[4] = nsave[6];
        printf("nddo:%12d", nsave[4]);
      }
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
      {
        nsave[5] = d->nddo;
        printf("nddo:%12d\n", nsave[5]);
      }
    } else {
#pragma omp single
      {
        nsave[5] = nsave[6];
        printf("nddo:%12d\n", nsave[5]);
      }
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
}

// void equilibrium_sc2(DEOMAUX *daux, DEOM *d, CTRL *c, const char *frho_name) {
//   vector<int> nsave(7);
//   FILE *frho = fopen(frho_name, "w");

//   const double dt = c->dt;
//   int nt_i, nt_f, nt;
//   double ti;
//   if (c->ti < 0) {
//     nt_i = ceil(abs(c->ti) / dt);
//     nt_f = ceil(abs(c->tf) / dt);
//     nt = nt_i + nt_f;
//     ti = -nt_i * dt;
//   } else {
//     nt_i = ceil(abs(c->ti) / dt);
//     nt_f = ceil(abs(c->tf) / dt);
//     nt = nt_f - nt_i;
//     ti = nt_i * dt;
//   }

//   int i;
//   int iado_;
//   const int n_threads = omp_get_max_threads();

//   printf("===   BRIEFING !!!!   ===\n");
//   printf("dt:%e\n", dt);
//   printf("nt:%d\n", nt);
//   printf("tf:%f\n", c->tf);
//   print(d);

//   printf("task: equilibrium, cal \\rho^{eq}\n");
//   printf("===   BRIEFING DONE   ===\n");

//   Trie *tree = new Trie(10 * d->nind);
//   tree->try_insert(1, 0);

//   const double dt2 = dt / 2.0;
//   const double dt3 = dt / 3.0;
//   const double dt4 = dt / 4.0;
//   const double dt6 = dt / 6.0;

//   d->hamt = d->ham1;
//   for (int mp = 0; mp < d->nind; mp++) {
//     int m1 = d->modLabel(mp);
//     d->qmdt[mp].noalias() = d->qmd1[m1];
//   }
//   int pos = tree->find(1);
// #pragma omp parallel default(shared)
//   for (int ii = 0; ii < nt; ii++) {
// #pragma omp single
//     {
//       fprintf(frho, "%16.6e", ti + ii * dt);
//       int pos = tree->find(1);
//       for (int i = 0; i < NSYS; i++)
//         fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());
//       fprintf(frho, "\n");

//       d->ddos[pos] = (d->ddos[pos] + d->ddos[pos].adjoint()) / 2;

//       printf("%16.6e\t", ti + ii * dt);
//       for (int i = 0; i < NSYS; i++)
//         printf("%16.6e\t", d->ddos[pos](i, i).real());
//       if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
//         printf("NH\t");
//       } else {
//         printf("HE\t");
//       }
//       // printf("%16.6e\t", d->ddos[pos](1, 0).real());
//       // printf("%16.6e\t", d->ddos[pos](1, 0).imag());

//       printf("nddo:%d\t", nsave[1]);
//       printf("nddo:%d\n", nsave[2]);
//     }

// #pragma omp for private(iado_) schedule(dynamic, 64)
//     for (iado_ = 0; iado_ < d->nddo; iado_++) {
//       for (int mp = 0; mp < 2 * d->nind + 1; mp++)
//         d->g_index_list_pos(iado_, mp) = 0;
//       for (int mp = 0; mp < d->nind; mp++)
//         daux->keys1(iado_, mp) = -1;
//       daux->ddos1[iado_].setZero();
//     }

// #pragma omp single
//     {
//       nsave[6] = d->nddo;
//       tree->clear();
//       delete tree;
//       tree = new Trie(min((ullint)d->nind * (ullint)d->nind * (ullint)d->nddo, d->nmax));
//       d->nddo = 0;
//       nsave[0] = 0;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = 0; i < nsave[6]; i++)
//       filter_p(d, d->ddos, d->keys, daux->ddos1, daux->keys1, tree, i);

// #pragma omp for private(iado_) schedule(dynamic, 64)
//     for (iado_ = 0; iado_ < nsave[6]; iado_++) {
//       d->keys.row(iado_) = daux->keys1.row(iado_);
//       d->ddos[iado_].noalias() = daux->ddos1[iado_];
//     }

// #pragma omp single
//     {
//       nsave[1] = d->nddo;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = nsave[0]; i < nsave[1]; i++)
//       construct_Mapping_p(d, tree, i);

// #pragma omp single
//     {
//       nsave[2] = d->nddo;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = 0; i < nsave[2]; i++)
//       syl_cal(daux->ddos1, d->ddos, d, i);

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = 0; i < nsave[2]; i++)
//       d->ddos[i].noalias() = daux->ddos1[i];
//   }
//   fclose(frho);
// }