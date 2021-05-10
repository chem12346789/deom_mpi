#include "deom.hpp"
#include "mr.hpp"

void equilibrium(DEOMAUX *daux, DEOM *d, CTRL *c, FILE *frho) {
  vector<int> nsave(7);
  MatrixXcd tmp;
  tmp = MatrixXcd::Zero(NSYS, NSYS);
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
  const int n_threads = omp_get_max_threads();

  printf("===   BRIEFING !!!!   ===\n");
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("n_threads:%d\n", n_threads);
  printf("nmax = %lld\n", d->nmax);
  printf("nsys = %d\n", d->nsys);
  printf("nind = %d\n", d->nind);
  printf("lmax = %d\n", d->lmax);
  printf("nham = %d\n", d->nham);
  printf("nmod = %d\n", d->nmod);
  printf("nmodmax = %d\n", d->nmodmax);
  printf("ferr = %g\n", d->ferr);
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

#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < d->nddo; iado_++) {
      for (int mp = 0; mp < 2 * d->nind + 1; mp++)
        daux->g_index_list_pos(iado_, mp) = 0;
      for (int mp = 0; mp < d->nind; mp++)
        daux->keys1(iado_, mp) = -1;
      daux->ddos1[iado_].setZero();
    }

#pragma omp single
    {
      nsave[6] = d->nddo;
      tree->clear();
      delete tree;
      tree = new Trie(min((ullint)2 * (ullint)d->nind * (ullint)d->nddo, d->nmax));
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
      construct_Mapping_p(d, d->keys, daux->g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[2] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[2]; i++)
      rem_cal(daux->ddos1, d->ddos, d, i, d->keys, daux->g_index_list_pos);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++)
      daux->ddos3[i].noalias() = d->ddos[i] + daux->ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[2]; i++)
      daux->ddos3[i].noalias() = daux->ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[2]; i++)
      construct_Mapping_p(daux->ddos3, d, d->keys, daux->g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[3] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[3]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i, d->keys, daux->g_index_list_pos);

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

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[3]; i++)
      construct_Mapping_p(daux->ddos3, d, d->keys, daux->g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[4] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[4]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i, d->keys, daux->g_index_list_pos);

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
    for (i = nsave[3]; i < nsave[4]; i++)
      for (int j = 0; j < NSYS; j++)
        for (int k = 0; k < NSYS; k++) {
          daux->ddos1[i].noalias() = daux->ddos2[i] * 2.0;
          daux->ddos3[i].noalias() = daux->ddos2[i] * dt;
        }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[4]; i++)
      construct_Mapping_p(daux->ddos3, d, d->keys, daux->g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[5] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[5]; i++)
      rem_cal(daux->ddos2, daux->ddos3, d, i, d->keys, daux->g_index_list_pos);

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
}

void equilibrium_sa(DEOMAUX *daux, DEOM *d, CTRL *c, FILE *frho) {
  DEOM *d1, deom1;
  CTRL *c1, ctrl1;
  d1 = &deom1;
  c1 = &ctrl1;

  equilibrium(daux, d, c, frho);
}