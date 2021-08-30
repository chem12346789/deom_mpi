#include "deom.hpp"
#include "index.hpp"

using namespace std;
using namespace json11;
using namespace Eigen;

void Init_dip(DIPOLE &dipole, DEOM *d, const Json &json) {
  dipole.sdip = nNNmat(d->nmodmax);
  for (int k = 0; k < d->nmodmax; k++)
    for (int i = 0; i < NSYS; i++)
      for (int j = 0; j < NSYS; j++)
        dipole.sdip[k](i, j) = json["sdip_cub"]["real"].array_items()[INDEX3(k, i, j)].number_value() + complex<double>(1i) * json["sdip_cub"]["imag"].array_items()[INDEX3(k, i, j)].number_value();

  dipole.pdip = nNNmat(d->nmodmax);
  for (int k = 0; k < d->nmodmax; k++)
    for (int i = 0; i < NSYS; i++)
      for (int j = 0; j < NSYS; j++)
        dipole.pdip[k](i, j) = json["pdip_cub"]["real"].array_items()[INDEX3(k, i, j)].number_value() + complex<double>(1i) * json["pdip_cub"]["imag"].array_items()[INDEX3(k, i, j)].number_value();

  dipole.bdip = VectorXcd::Zero(d->nind);
  for (int i = 0; i < d->nind; i++)
    dipole.bdip(i) = json["bdip_cub"]["real"].array_items()[i].number_value() + complex<double>(1i) * json["bdip_cub"]["imag"].array_items()[i].number_value();
  dipole.on = true;
}

void Init_pulse(PULSE &pulse, const Json &json) {
  // ## init pulse##
  pulse.on = true;
  pulse.ampl = json["ampl"].number_value();
  pulse.sigm = json["sigm"].number_value();
  pulse.freq = json["freq"].number_value();
  pulse.varp = json["varp"].number_value();
  // ## init pulse##
}

void Init_sys(DEOM *d, const Json &json, const char *input_name) {
  // ##  init para  ##
  d->nmax = json["nmax"].int_value();
  d->nsys = json["nsys"].int_value();
  d->nind = json["nind"].int_value();
  d->lmax = json["lmax"].int_value();
  d->nham = json["nham"].int_value();
  d->nmod = json["nmod"].int_value();
  d->nmodmax = json["nmodmax"].int_value();

  if (json["ferr"].is_number()) {
    d->ferr = json["ferr"].number_value();
  } else {
    d->ferr = 1e-14;
  }

  if (json["load"].is_number()) {
    d->load = json["load"].int_value();
  } else {
    d->load = 2;
  }

  printf("load coef: %d\n", d->load);
  // ##  init comb   ##
  d->combmax = d->nind + d->lmax + 2;
  d->comb_list.resize(d->combmax, d->combmax);
  for (int j = 0; j < d->combmax; j++)
    d->comb_list(0, j) = 0;
  d->comb_list(0, 0) = 1;

  for (int i = 1; i < d->combmax; i++) {
    for (int j = 1; j < d->combmax; j++) {
      d->comb_list(i, j) = d->comb_list((i - 1), j) + d->comb_list((i - 1), j - 1);
    }
    d->comb_list(i, 0) = 1;
  }
  // ##  init comb done  ##

  if (d->comb_list((d->lmax + d->nind), d->lmax) < d->nmax) {
    d->nmax = d->comb_list(d->lmax + d->nind, d->lmax);
  } else {
    printf("mem error! but do not panic it maybe right.\n");
  }
  // ##  init para done ##

  if (NSYS != d->nsys) {
    printf("ERROR NYSY!! PEALASE REPLACE NSYS AND RECOMPILE IT!!!\nBTW: the original design is we do not change it during the workflow, and fixed mat will gain about 2.0x speed up, so we fix it, if you prefer another feature just change it.");
    exit(2);
  }

  if (d->lmax > 10000) {
    printf("It will cost too much, consider change another system or method!!");
    exit(3);
  }

  // ##  init list things ##
  d->ham1 = MatrixXcd::Zero(NSYS, NSYS);
  for (int i = 0; i < d->nham; i++)
    for (int j = 0; j < d->nham; j++)
      d->ham1(i, j) = json["ham1"]["real"].array_items()[INDEX2(i, j, d->nham)].number_value() + complex<double>(1i) * json["ham1"]["imag"].array_items()[INDEX2(i, j, d->nham)].number_value();
  std::cout << d->ham1 << std::endl;

  d->modLabel = VectorXi(d->nind);
  for (int i = 0; i < d->nind; i++)
    d->modLabel(i) = json["modLabel"]["real"].array_items()[i].int_value();

  d->delt_res = VectorXcd::Zero(d->nmodmax);

  d->coef_abs = VectorXcd::Zero(d->nind);
  d->expn = VectorXcd::Zero(d->nind);
  for (int i = 0; i < d->nmodmax; i++)
    d->delt_res(i) = json["delt_res"]["real"].array_items()[i].number_value() + complex<double>(1i) * json["delt_res"]["imag"].array_items()[i].number_value();

  for (int i = 0; i < d->nmodmax; i++)
    if (abs(d->delt_res(i)) > 1e-15) {
      printf("PLZ use pade = 1!!!!\nIt can save u a lot of time.");
      exit(201);
    }

  for (int i = 0; i < d->nind; i++) {
    d->coef_abs(i) = sqrt(json["coef_abs"]["real"].array_items()[i].number_value() + complex<double>(1i) * json["coef_abs"]["imag"].array_items()[i].number_value());
    d->expn(i) = json["expn"]["real"].array_items()[i].number_value() + complex<double>(1i) * json["expn"]["imag"].array_items()[i].number_value();
  }

  d->coef_lft = std::vector<MatrixXcd, Eigen::aligned_allocator<MatrixXcd>>(d->nind);
  d->coef_rht = std::vector<MatrixXcd, Eigen::aligned_allocator<MatrixXcd>>(d->nind);
  for (int mp = 0; mp < d->nind; mp++) {
    d->coef_lft[mp].resize(d->nmodmax, d->nmodmax);
    d->coef_rht[mp].resize(d->nmodmax, d->nmodmax);
    for (int ii = 0; ii < d->nmodmax; ii++) {
      for (int m = 0; m < d->nmodmax; m++) {
        d->coef_lft[mp](m, ii) = json["coef_lft"]["real"].array_items()[INDEX3(mp, m, ii, d->nmodmax)].number_value() + complex<double>(1i) * json["coef_lft"]["imag"].array_items()[INDEX3(mp, m, ii, d->nmodmax)].number_value();
        d->coef_rht[mp](m, ii) = json["coef_rht"]["real"].array_items()[INDEX3(mp, m, ii, d->nmodmax)].number_value() + complex<double>(1i) * json["coef_rht"]["imag"].array_items()[INDEX3(mp, m, ii, d->nmodmax)].number_value();
      }
    }
  }

  d->qmd1 = nNNmat(d->nmodmax);
  for (int i = 0; i < d->nmodmax; i++)
    for (int j = 0; j < NSYS; j++)
      for (int k = 0; k < NSYS; k++)
        d->qmd1[i](j, k) = json["qmd1"]["real"].array_items()[INDEX3(i, j, k)].number_value() + complex<double>(1i) * json["qmd1"]["imag"].array_items()[INDEX3(i, j, k)].number_value();
  // ##  init list/array things done ##

  d->g_index_list_pos = XXimat_r(d->nmax, (2 * d->nind + 1));
  for (int iado = 0; iado < d->nmax; iado++)
    for (int mp = 0; mp < 2 * d->nind + 1; mp++)
      d->g_index_list_pos(iado, mp) = 0;

  // ##  init keys and gams ##
  d->keys = XXuimat_r(d->nmax, d->nind);
  for (int iado = 0; iado < d->nmax; iado++) {
    for (int mp = 0; mp < d->nind; mp++) {
      d->keys(iado, mp) = 0;
    }
  }

  d->zerokey = XXuivet_r(d->nind);
  for (int i = 0; i < d->nind; i++) {
    d->zerokey[i] = 0;
  }

  d->keys.row(0) = d->zerokey;
  // ##  init qmdt and hamt  ##
  d->hamt = MatrixNcd::Zero();
  d->qmdt = nNNmat(d->nind);
  // ##  init qmdt and hamt done ##

  d->ddos = nNNmat(d->nmax);
  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    d->ddos[iado_].setZero();
  }

  d->nddo = 1;
  d->ddos[0](0, 0) = 1;
}

void Init_ctrl(CTRL *c, const Json &json) {
  c->dt = json["dt"].number_value();
  c->ti = json["ti"].number_value();
  c->tf = json["tf"].number_value();
}

void Init_ctrl(CTRL *c, double dt, double ti, double tf) {
  c->dt = dt;
  c->ti = ti;
  c->tf = tf;
}

void Init_aux(DEOMAUX *daux, DEOM *d, const Json &json) {
  daux->keys1 = XXuimat_r(d->nmax, d->nind);

  daux->key_tmp = XXuivet_r(d->nind);
  // // daux->key_tmp1 = VectorXi(d->nind);
  // // for (int key_index = 0; key_index < d->nind; key_index++)
  // //   daux->key_tmp(key_index) = 0;

  // ##  init keys and gams done  ##
  daux->ddos1 = nNNmat(d->nmax);
  daux->ddos2 = nNNmat(d->nmax);
  daux->ddos3 = nNNmat(d->nmax);
}