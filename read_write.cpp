#include "deom.hpp"
#include "mr.hpp"

void write_ddos(DEOMAUX *daux, DEOM *d, const char *fname_ddos, const char *fname_keys) {
  ofstream fstream_ddos(fname_ddos);
  ofstream fstream_keys(fname_keys);
  IOFormat HeavyFmt(FullPrecision);

  int i;
  int iado_;
  vector<int> nsave(7);
  Trie *tree = new Trie(10 * d->nind);

#pragma omp parallel default(shared)
  {
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
  }
  for (size_t i = 0; i < d->nddo; i++) {
    fstream_ddos << d->ddos[i].format(HeavyFmt) << '\n';
  }
  fstream_keys << d->keys.block(0, 0, d->nddo, d->nind);
  fstream_ddos.close();
  fstream_keys.close();
}