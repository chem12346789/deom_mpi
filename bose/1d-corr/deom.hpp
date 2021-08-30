#pragma once
#include "headfile.hpp"

typedef struct {
  bool on;
  nNNmat sdip;
  nNNmat pdip;
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
  int load;
  int nmod;
  int nmodmax;
  int nham;
  VectorXi modLabel;
  MatrixNcd ham1;
  MatrixNcd hamt;
  nNNmat qmd1;
  VectorXcd delt_res;
  VectorXcd expn;
  VectorXcd coef_abs;
  nXXmat coef_lft;
  nXXmat coef_rht;
  XXuimat_r keys;
  DIPOLE dipole;
  DIPOLE dipole1;
  nullmat comb_list;
  nNNmat qmdt;
  XXuivet_r zerokey;
  nNNmat ddos;
  PULSE pulse;
  XXimat_r g_index_list_pos;
} DEOM;

typedef struct {
  double dt;
  double ti;
  double tf;
} CTRL;

typedef struct {
  XXuivet_r key_tmp;
  XXuivet_r key_tmp1;
  XXuimat_r keys1;
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
  printf("##  ham1\n");
  std::cout << fixed << showpos << d->ham1.format(HeavyFmt) << std::endl;
  printf("##  ham1 done\n\n");

  printf("##  qmd1\n");
  for (int i = 0; i < d->nmodmax; i++)
    std::cout << fixed << showpos << d->qmd1[i].format(HeavyFmt) << std::endl;
  printf("##  qmd1 done\n\n");

  printf("##  expn\n");
  std::cout << fixed << showpos << d->expn.format(HeavyFmt) << std::endl;
  printf("##  expn done\n\n");

  printf("##  coef_abs\n");
  std::cout << fixed << showpos << d->coef_abs.format(HeavyFmt) << std::endl;
  printf("##  coef_abs done\n\n");

  printf("##  coef_lft\n");
  for (int mp = 0; mp < d->nind; mp++)
    std::cout << fixed << showpos << d->coef_lft[mp].format(HeavyFmt) << std::endl;
  printf("##  coef_lft done\n\n");

  printf("##  coef_rht\n");
  for (int mp = 0; mp < d->nind; mp++)
    std::cout << fixed << showpos << d->coef_rht[mp].format(HeavyFmt) << std::endl;
  printf("##  coef_rht done\n\n");

  printf("##  delt_res\n");
  std::cout << fixed << showpos << d->delt_res.format(HeavyFmt) << std::endl;
  printf("##  delt_res done\n\n");
}

void print(DIPOLE dipole, DEOM *d) {
  IOFormat HeavyFmt(StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

  printf("##  sdip\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << fixed << showpos << d->dipole.sdip[k].format(HeavyFmt) << std::endl;
  }
  printf("##  sdip done\n\n");
  printf("##  pdip\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << fixed << showpos << d->dipole.pdip[k].format(HeavyFmt) << std::endl;
  }
  printf("##  pdip done\n\n");
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
