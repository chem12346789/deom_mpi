#include "deom.hpp"
#include "mr.hpp"

void strong_field(DEOMAUX *daux, DEOM *d, CTRL *c, FILE *frho, FILE *fpolar) {
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

  print(d);
  if (d->pulse.on) {
    print(d->pulse);
  }

  if (d->dipole.on) {
    print(d->dipole, d);
  }


  printf("task: cal strong_field <\\P(t)\\P(0)>\n");
  printf("===   BRIEFING DONE   ===\n");

  Trie *tree = new Trie(10 * d->nind);
  tree->try_insert(1, 0);
  const double dt2 = dt / 2.0;
  const double dt3 = dt / 3.0;
  const double dt4 = dt / 4.0;
  const double dt6 = dt / 6.0;

#pragma omp parallel default(shared)
  for (int ii = 0; ii < nt; ii++) {
#pragma omp single
    {
      fprintf(frho, "%16.6e\n", ti + ii * dt);
      int pos = tree->find(1);
      tmp.setZero();
      for (int nmod = 0; nmod < d->nmodmax; nmod++)
        tmp += d->dipole.sdip[nmod] * d->ddos[pos];

      for (int mp = 0; mp < d->nind; mp++) {
        const complex<double> sn = d->dipole.bdip(mp) * d->coef_abs(mp);
        gen_key(daux->key_tmp, daux->key_tmp1, mp, 1, d);
        int m1 = d->modLabel(mp);
        pos = tree->find(hash_bad(daux->key_tmp1, d));
        if (pos > -1) {
          tmp += sn * d->dipole.pdip[m1] * d->ddos[pos];
        }
      }

      complex<double> trace = tmp.trace();

      fprintf(fpolar, "%16.6e", ti + ii * dt);
      fprintf(fpolar, "%20.10e%20.10e", real(trace), imag(trace));
      fprintf(fpolar, "\n");
      pos = tree->find(1);

      for (int i = 0; i < NSYS; i++)
        fprintf(frho, "%16.6e\t", d->ddos[pos](i, i).real());

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
      construct_Mapping_p(d, d->keys, daux->g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[2] = d->nddo;
      double et = pulse_et(ti + ii * dt, d->pulse);

      d->hamt = d->ham1;
      for (int nmod = 0; nmod < d->nmodmax; nmod++)
        (d->hamt).noalias() -= d->dipole.sdip[nmod] * et;

      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        d->qmdt[mp].noalias() = d->qmd1[m1] - d->dipole.pdip[m1] * d->dipole.bdip(mp) * et;
      }
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
      double et = pulse_et(ti + ii * dt + dt2, d->pulse);
      d->hamt = d->ham1;
      for (int nmod = 0; nmod < d->nmodmax; nmod++)
        (d->hamt).noalias() -= d->dipole.sdip[nmod] * et;

      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        d->qmdt[mp].noalias() = d->qmd1[m1] - d->dipole.pdip[m1] * d->dipole.bdip(mp) * et;
      }
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
      double et = pulse_et(ti + ii * dt + dt2, d->pulse);
      d->hamt = d->ham1;
      for (int nmod = 0; nmod < d->nmodmax; nmod++)
        (d->hamt).noalias() -= d->dipole.sdip[nmod] * et;

      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        d->qmdt[mp].noalias() = d->qmd1[m1] - d->dipole.pdip[m1] * d->dipole.bdip(mp) * et;
      }
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
      double et = pulse_et(ti + ii * dt + dt, d->pulse);
      d->hamt = d->ham1;
      for (int nmod = 0; nmod < d->nmodmax; nmod++)
        (d->hamt).noalias() -= d->dipole.sdip[nmod] * et;

      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        d->qmdt[mp].noalias() = d->qmd1[m1] - d->dipole.pdip[m1] * d->dipole.bdip(mp) * et;
      }
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
