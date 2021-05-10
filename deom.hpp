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
  int nmod;
  int nmodmax;
  int nham;
  VectorXi modLabel;
  Matrix2cd ham1;
  VectorXcd delt_res;
  nNNmat qmd1;
  XXimat_r keys;
  VectorXcd coef_abs;
  nXXmat coef_lft;
  nXXmat coef_rht;
  VectorXcd expn;
  DIPOLE dipole;
  nullmat comb_list;
  nNNmat qmdt;
  Matrix2cd hamt;
  VectorXi emptykey;
  nNNmat ddos;
  PULSE pulse;
} DEOM;

typedef struct {
  double dt;
  double ti;
  double tf;
} CTRL;

typedef struct {
  VectorXi key_tmp;
  VectorXi key_tmp1;
  XXimat_r g_index_list_pos;
  XXimat_r keys1;
  nNNmat ddos1;
  nNNmat ddos2;
  nNNmat ddos3;
} DEOMAUX;

void print(DEOM *d) {
  
}

void print(DIPOLE dipole, DEOM *d) {
  printf("##  sdip\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << d->dipole.sdip[k] << std::endl;
  }
  printf("##  sdip done\n\n");
  printf("##  pdip\n");
  for (int k = 0; k < d->nmodmax; k++) {
    printf("slice %d\n", k);
    std::cout << d->dipole.pdip[k] << std::endl;
  }
  printf("##  pdip done\n\n");
  printf("##  bdip\n");
  std::cout << d->dipole.bdip << std::endl;
  printf("##  bdip done\n\n");
}

void print(PULSE pulse) {
  printf("ampl\t%e\n", pulse.ampl);
  printf("sigm\t%e\n", pulse.sigm);
  printf("freq\t%e\n", pulse.freq);
  printf("varp\t%e\n", pulse.varp);
}