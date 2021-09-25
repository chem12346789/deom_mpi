#include "deom.hpp"
#include "mr.hpp"
// Calculate equilibrium state. You will find in small system, equilibrium_sc2 is faster than equilibrium.
// The target is equilibrium state, but you can use it to do some "propgator" in \rho.
// i/o parameter:
// frho_name: each step of ddos0
// some contorl parameter:
// ti: init time
// tf: fimal time
// dt: delta time
// re_tree: step to rebuild the tree
// filter: step to do filter
void equilibrium(DEOMAUX *daux, DEOM *d, CTRL *c, const char *frho_name) {
  vector<int> nsave(7);
  FILE *frho = fopen(frho_name, "w");

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

  printf("\t===   BRIEFING !!!!   ===\n");
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("tf:%f\n", c->tf);
  print(d);

  printf("task: equilibrium, fokker_planck, not standard deom, cal \\rho^{eq}\n");
  printf("\t===   BRIEFING DONE   ===\n");

  Trie *tree = new Trie(10 * d->nind);
  const double dt2 = dt / 2.0;
  const double dt3 = dt / 3.0;
  const double dt4 = dt / 4.0;
  const double dt6 = dt / 6.0;
  int re_tree = c->re_tree==0?1:c->re_tree;
  int filter = c->filter==0?1:c->filter;
  ullint tree_size;

  d->hamt = d->ham1;
  for (int mp = 0; mp < d->nind; mp++) {
    int m1 = d->modLabel(mp);
    d->qmdt[mp].noalias() = d->qmd1[m1];
  }

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

      if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
        printf("NH\n");
      } else {
        printf("HE\n");
      }

      nsave[6] = d->nddo;
      ullint tree_size_n = tree->size();
      if (ii < re_tree) {
        tree->clear();
        tree_size = min(5 * (ullint)d->nind * (ullint)d->nddo, d->nmax);
        tree = new Trie(tree_size);
      } else {
        if (ii % re_tree == 0) {
          tree->clear();
          delete tree;
          tree_size = min((ullint)d->nddo + (ullint)(0.5 * d->nddo), (ullint)2 * d->nmax);
          tree = new Trie(tree_size);
        }
      }

      if (ii < filter || ii % filter == 0) {
        printf("set tree size: %llu\t real tree size: %llu\t", tree_size, tree_size_n);
        d->nddo = 0;
        nsave[0] = 0;
      }
    }

    if (ii < filter || ii % filter == 0) {
#pragma omp for private(iado_) schedule(dynamic, 64)
      for (iado_ = 0; iado_ < nsave[6]; iado_++) {
        d->g_index_list_pos.row(iado_).setZero();
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
        printf("nddo:%12d", nsave[1]);
        printf("nddo:%12d", nsave[2]);
      }
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

// void equilibrium(DEOMAUX *daux, DEOM *d, CTRL *c, const char *frho_name) {
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
//   ullint tree_size;

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

// #pragma omp parallel default(shared)
//   for (int ii = 0; ii < nt; ii++) {
// #pragma omp single
//     {
//       fprintf(frho, "%16.6e", ti + ii * dt);
//       int pos = tree->find(1);
//       for (int i = 0; i < NSYS; i++)
//         fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());
//       fprintf(frho, "\n");

//       printf("%16.6e\t", ti + ii * dt);
//       for (int i = 0; i < NSYS; i++)
//         printf("%16.6e\t", d->ddos[pos](i, i).real());

//       if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
//         printf("NH\t");
//       } else {
//         printf("HE\t");
//       }

//       printf("nddo:%d\t", nsave[1]);
//       printf("nddo:%d\t", nsave[2]);
//       printf("nddo:%d\t", nsave[3]);
//       printf("nddo:%d\t", nsave[4]);
//       printf("nddo:%d\n", nsave[5]);
//     }

// #pragma omp for private(iado_) schedule(dynamic, 64)
//     for (iado_ = 0; iado_ < d->nddo; iado_++) {
//       d->g_index_list_pos.row(iado_).setZero();
//       // daux->ddos1[iado_].setZero();
//     }

// #pragma omp single
//     {
//       nsave[6] = d->nddo;
//       ullint tree_size_n = tree->size();
//       tree->clear();
//       if (ii < 50) {
//         tree_size = min((ullint)d->nind * (ullint)d->nind * (ullint)d->nddo, d->nmax);
//         tree = new Trie(tree_size);
//       } else {
//         if (ii % 50 == 0) {
//           delete tree;
//           tree_size = min((ullint)d->nddo + (ullint)(1.5 * d->nddo), (ullint)2 * d->nmax);
//           tree = new Trie(tree_size);
//         }
//       }

//       printf("set tree size: %llu\t real tree size: %llu\t", tree_size, tree_size_n);
//       d->nddo = 0;
//       nsave[0] = 0;
//     }
//     // {
//     //   nsave[6] = d->nddo;
//     //   tree->clear();
//     //   delete tree;
//     //   tree = new Trie(min((ullint)d->load * (ullint)d->nind * (ullint)d->nddo, d->nmax));
//     //   d->nddo = 0;
//     //   nsave[0] = 0;
//     // }

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
//       rem_cal(daux->ddos1, d->ddos, d, i);

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = 0; i < nsave[1]; i++)
//       daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos1[i] * dt2;

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[1]; i < nsave[2]; i++)
//       daux->ddos3[i].noalias() = daux->ddos1[i] * dt2;

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = nsave[0]; i < nsave[2]; i++)
//       construct_Mapping_p(daux->ddos3, d, tree, i);

// #pragma omp single
//     {
//       nsave[3] = d->nddo;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = 0; i < nsave[3]; i++)
//       rem_cal(daux->ddos2, daux->ddos3, d, i);

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = 0; i < nsave[1]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos2[i] * dt2;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[1]; i < nsave[2]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = daux->ddos2[i] * dt2;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[2]; i < nsave[3]; i++) {
//       daux->ddos1[i].noalias() = daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = daux->ddos2[i] * dt2;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = nsave[0]; i < nsave[3]; i++)
//       construct_Mapping_p(daux->ddos3, d, tree, i);

// #pragma omp single
//     {
//       nsave[4] = d->nddo;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = 0; i < nsave[4]; i++)
//       rem_cal(daux->ddos2, daux->ddos3, d, i);

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = 0; i < nsave[1]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos2[i] * dt;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[1]; i < nsave[3]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = daux->ddos2[i] * dt;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[3]; i < nsave[4]; i++) {
//       daux->ddos1[i].noalias() = daux->ddos2[i] * 2.0;
//       daux->ddos3[i].noalias() = daux->ddos2[i] * dt;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = nsave[0]; i < nsave[4]; i++)
//       construct_Mapping_p(daux->ddos3, d, tree, i);

// #pragma omp single
//     {
//       nsave[5] = d->nddo;
//     }

// #pragma omp for private(i) schedule(dynamic, 16)
//     for (i = 0; i < nsave[5]; i++)
//       rem_cal(daux->ddos2, daux->ddos3, d, i);

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = 0; i < nsave[1]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i];
//       d->ddos[i].noalias() += daux->ddos1[i] * dt6;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[1]; i < nsave[4]; i++) {
//       daux->ddos1[i].noalias() += daux->ddos2[i];
//       d->ddos[i].noalias() = daux->ddos1[i] * dt6;
//     }

// #pragma omp for private(i) schedule(dynamic, 64)
//     for (i = nsave[4]; i < nsave[5]; i++) {
//       daux->ddos1[i].noalias() = daux->ddos2[i];
//       d->ddos[i].noalias() = daux->ddos1[i] * dt6;
//     }
//   }

//   fclose(frho);
// }

void equilibrium_sc2(DEOMAUX *daux, DEOM *d, CTRL *c, const char *frho_name) {
  vector<int> nsave(7);
  FILE *frho = fopen(frho_name, "w");

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

  printf("===   BRIEFING !!!!   ===\n");
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("tf:%f\n", c->tf);
  print(d);

  printf("task: equilibrium, cal \\rho^{eq}\n");
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
    d->qmdt[mp].noalias() = d->qmd1[m1];
  }
  int pos = tree->find(1);
#pragma omp parallel default(shared)
  for (int ii = 0; ii < nt; ii++) {
#pragma omp single
    {
      fprintf(frho, "%16.6e", ti + ii * dt);
      int pos = tree->find(1);
      for (int i = 0; i < NSYS; i++)
        fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());
      fprintf(frho, "\n");

      d->ddos[pos] = (d->ddos[pos] + d->ddos[pos].adjoint()) / 2;

      printf("%16.6e\t", ti + ii * dt);
      for (int i = 0; i < NSYS; i++)
        printf("%16.6e\t", d->ddos[pos](i, i).real());
      if (is_valid(d->ddos[pos] - d->ddos[pos].adjoint())) {
        printf("NH\t");
      } else {
        printf("HE\t");
      }
      // printf("%16.6e\t", d->ddos[pos](1, 0).real());
      // printf("%16.6e\t", d->ddos[pos](1, 0).imag());

      printf("nddo:%d\t", nsave[1]);
      printf("nddo:%d\n", nsave[2]);
    }

#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < d->nddo; iado_++) {
      d->g_index_list_pos.row(iado_).setZero();
      // daux->ddos1[iado_].setZero();
    }

#pragma omp single
    {
      nsave[6] = d->nddo;
      tree->clear();
      delete tree;
      tree = new Trie(min((ullint)d->load * (ullint)d->nind * (ullint)d->nddo, d->nmax));
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
      syl_cal(daux->ddos1, d->ddos, d, i);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[2]; i++)
      d->ddos[i].noalias() = daux->ddos1[i];
  }
  fclose(frho);
}