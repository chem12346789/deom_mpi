#include "deom.hpp"
#include "mr.hpp"

// this file is same at both fermion and boson keep it in mind

// write ddos and keys to file
// we will write ddos to file(and read it later)
// this function should NOT change any ddos or keys in d(daux may be changed)
void write_ddos(DEOMAUX *daux, DEOM *d, const char *fname_ddos, const char *fname_keys) {
  ofstream fstream_ddos(fname_ddos);
  ofstream fstream_keys(fname_keys);
  IOFormat HeavyFmt(FullPrecision);

  for (size_t i = 0; i < d->nddo; i++) {
    if(is_valid(d->ddos[i], d->ferr)){
    fstream_ddos << d->ddos[i].format(HeavyFmt) << '\n';
    fstream_keys << d->keys.row(d->nddo);}
  }
  fstream_ddos.close();
  fstream_keys.close();
}

// write ddos and keys to file
// we will write ddos to three part file(real one , imag one and both of it)
// so that you can save it for python read later
void write_ddos(DEOMAUX *daux, DEOM *d, const char *fname_ddos, const char *fname_ddosr, const char *fname_ddosi, const char *fname_keys) {
  ofstream fstream_ddosr(fname_ddosr);
  ofstream fstream_ddosi(fname_ddosi);
  ofstream fstream_ddos(fname_ddos);
  ofstream fstream_keys(fname_keys);
  IOFormat HeavyFmt(FullPrecision);

  int i;
  int iado_;
  vector<int> nsave(7);
  Trie *tree;

#pragma omp parallel default(shared)
  {
#pragma omp single
    {
      nsave[6] = d->nddo;
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
  }

  for (size_t i = 0; i < d->nddo; i++) {
    fstream_ddosr << d->ddos[i].real().format(HeavyFmt) << '\n';
    fstream_ddosi << d->ddos[i].imag().format(HeavyFmt) << '\n';
    fstream_ddos << d->ddos[i].format(HeavyFmt) << '\n';
  }
  fstream_keys << d->keys.block(0, 0, d->nddo, d->nind);
  fstream_ddosr.close();
  fstream_ddosi.close();
  fstream_ddos.close();
  fstream_keys.close();
}

// read ddos, you can load the ddos file you save before. NOTE you should init deom before doing that.
void read_ddos(DEOMAUX *daux, DEOM *d, const char *fname_ddos, const char *fname_keys) {
  ifstream fstream_ddos(fname_ddos);
  ifstream fstream_keys(fname_keys);

  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    d->keys.row(iado_).setZero();
    d->ddos[iado_].setZero();
  }
  
  int nrow_keys = 0;
  while (!fstream_keys.eof()) {
    string line;
    getline(fstream_keys, line);
    stringstream stream(line);
    for (int ncol = 0; ncol < d->nind; ncol++)
      stream >> d->keys(nrow_keys, ncol);
    nrow_keys += 1;
  }

  int iado = 0;
  while (!fstream_ddos.eof()) {
    string line;
    for (int nrow = 0; nrow < d->nsys; nrow++) {
      getline(fstream_ddos, line);
      stringstream stream(line);
      for (int ncol = 0; ncol < d->nsys; ncol++)
        stream >> d->ddos[iado](nrow, ncol);
    }
    if (!(line.empty()))
      iado += 1;
  }

  printf("%d\n", nrow_keys);
  printf("%d\n", iado);
  if (iado == nrow_keys)
    d->nddo = nrow_keys;
  else
    exit(202);
  fstream_ddos.close();
  fstream_keys.close();
}

