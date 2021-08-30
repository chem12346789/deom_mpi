#pragma once
#include "headfile.hpp"

typedef struct {
  bool on;
  nNNmat sdip;
  nNNmat pdipa;
  nNNmat pdipc;
  VectorXcd bdip;
} DIPOLE;

typedef struct {
  ullint nmax;
  int nsys;
  int nddo;
  float ferr;
  int combmax;
  int nind;
  int lmax;
  int nmod;
  int nmodmax;
  int nham;
  VectorXi modLabel;
  VectorXi nmodmaxLabel;
  MatrixNcd ham1;
  MatrixNcd hamt;
  nNNmat qmd1a;
  nNNmat qmd1c;
  VectorXcd delt_res;
  VectorXcd expn;
  VectorXcd coef_abs;
  nXXmat coef_lft;
  nXXmat coef_rht;
  XXbmat_r keys;
  DIPOLE dipole;
  DIPOLE dipole1;
  nullmat comb_list;
  nNNmat qmdta;
  nNNmat qmdtc;
  XXbvet_r zerokey;
  XXbvet_r emptykey;
  nNNmat ddos;
  PULSE pulse;
  string input_name;
  string input_hash;
  XXimat_r g_index_list_pos;
} DEOM;

typedef struct {
  double dt;
  double ti;
  double tf;
} CTRL;

typedef struct {
  XXbvet_r key_tmp;
  XXbvet_r key_tmp1;
  XXbmat_r keys1;
  nNNmat ddos1;
  nNNmat ddos2;
  nNNmat ddos3;
} DEOMAUX;

void print(DEOM *d) {
  IOFormat HeavyFmt(StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

  printf("nmax = %lld\n", d->nmax);
  printf("nsys = %d\n", d->nsys);
  printf("nind = %d\n", d->nind);
  printf("lmax = %d\n", d->lmax);
  printf("nham = %d\n", d->nham);
  printf("nmod = %d\n", d->nmod);
  printf("nmodmax = %d\n", d->nmodmax);
  printf("ferr = %g\n", d->ferr);
  printf("##  system hamiltonian\n");
  std::cout << fixed << showpos << d->ham1.format(HeavyFmt) << std::endl;
  printf("##  system hamiltonian done\n\n");

  printf("##  creation and annihilation operators(sigma first)\n");
  for (int i = 0; i < d->nmod; i++)
    std::cout << fixed << showpos << d->qmd1a[i].format(HeavyFmt) << std::endl;
  printf("##  creation and annihilation operators(sigma first) done\n\n");

  printf("##  creation and annihilation operators(sigma bar first)\n");
  for (int i = 0; i < d->nmod; i++)
    std::cout << fixed << showpos << d->qmd1c[i].format(HeavyFmt) << std::endl;
  printf("##  creation and annihilation operators(sigma bar first) done\n\n");

  printf("##  lamba in correlation function\n");
  std::cout << fixed << showpos << d->expn.format(HeavyFmt) << std::endl;
  printf("##  lamba in correlation function done\n\n");

  printf("##  modify coefficient eta in correlation function\n");
  std::cout << fixed << showpos << d->coef_abs.format(HeavyFmt) << std::endl;
  printf("##  modify coefficient eta in correlation function done\n\n");

  printf("##  eta in correlation function\n");
  for (int mp = 0; mp < d->nind; mp++)
    std::cout << fixed << showpos << d->coef_lft[mp].format(HeavyFmt) << std::endl;
  printf("##  eta in correlation function done\n\n");

  printf("##  eta star in correlation function\n");
  for (int mp = 0; mp < d->nind; mp++)
    std::cout << fixed << showpos << d->coef_rht[mp].format(HeavyFmt) << std::endl;
  printf("##  eta star in correlation function done\n\n");

  printf("##  excessive pade note it will be zero\n");
  std::cout << fixed << showpos << d->delt_res.format(HeavyFmt) << std::endl;
  printf("##  excessive pade done\n\n");

  printf("##  mod label to distinguish different pade or drude mode\n");
  std::cout << fixed << showpos << d->modLabel.format(HeavyFmt) << std::endl;
  printf("##  mod label done\n\n");

  printf("##  mod label to distinguish different mode in the same pade or drude mode\n");
  std::cout << fixed << showpos << d->modLabel.format(HeavyFmt) << std::endl;
  printf("##  mod label done\n\n");
}

void print(DIPOLE dipole, DEOM *d) {
  IOFormat HeavyFmt(StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

  printf("##  sdip\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << fixed << showpos << d->dipole.sdip[k].format(HeavyFmt) << std::endl;
  }
  printf("##  sdip done\n\n");
  printf("##  pdipa\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << fixed << showpos << d->dipole.pdipa[k].format(HeavyFmt) << std::endl;
  }
  printf("##  pdipa done\n\n");

  printf("##  pdipc\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << fixed << showpos << d->dipole.pdipc[k].format(HeavyFmt) << std::endl;
  }
  printf("##  pdipc done\n\n");

  printf("##  bdip\n");
  std::cout << fixed << showpos << d->dipole.bdip.format(HeavyFmt) << std::endl;
  printf("##  bdip done\n\n");
}

void print(PULSE pulse) {
  printf("ampl\t%e\n", pulse.ampl);
  printf("sigm\t%e\n", pulse.sigm);
  printf("freq\t%e\n", pulse.freq);
  printf("varp\t%e\n", pulse.varp);
}
