// #include "algebra.cpp"
#include "deom.hpp"
#include "mr.cpp"
// #include "pulse.hpp"
#include <gperftools/profiler.h>
#include <stdio.h>
#include <thread> // std::thread
#include <unistd.h>
using namespace std;
using Eigen::Matrix2cd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void init_deom(const Json &json, int argc, char *argv[]) {
  printf("init\n");
  DEOM *d, deom;
  d = &deom;

  // ##  init para  ##
  d->nmax = json["nmax"].int_value();
  d->nsys = json["nsys"].int_value();
  d->nind = json["nind"].int_value();
  d->lmax = json["lmax"].int_value();
  d->nham = json["nham"].int_value();
  d->nmod = json["nmod"].int_value();
  d->nmodmax = json["nmodmax"].int_value();
  d->ferr = json["ferr"].number_value();

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
  printf("nmax = %lld\n", d->nmax);
  printf("nsys = %d\n", d->nsys);
  printf("nind = %d\n", d->nind);
  printf("lmax = %d\n", d->lmax);
  printf("nham = %d\n", d->nham);
  printf("nmod = %d\n", d->nmod);
  printf("nmodmax = %d\n", d->nmodmax);
  printf("ferr = %g\n", d->ferr);
  // ##  init para done ##

  if (NSYS != d->nsys) {
    printf("ERROR NYSY!! PEALASE REPLACE NSYS AND RECOMPILE IT!!!\nBTW: the original design is we do not change it during the workflow, and fixed mat will gain about 2.0x speed up, so we fix it, if you prefer another feature just change it.");
    exit(1);
  }

  if (d->nmodmax != 1) {
    printf("ERROR nmodmax, this code is NOT contain nul mode!!!");
    exit(2);
  }

  // ##  init list things ##
  d->ham1 = MatrixXcd::Zero(NSYS, NSYS);
  for (int i = 0; i < d->nham; i++)
    for (int j = 0; j < d->nham; j++)
      d->ham1(i, j) = json["ham1"]["real"].array_items()[INDEX2(i, j, d->nham)].number_value() + complex<double>(1i) * json["ham1"]["imag"].array_items()[INDEX2(i, j, d->nham)].number_value();
  printf("# ham1\n");
  std::cout << d->ham1 << std::endl;

  d->modLabel = VectorXi(d->nind);
  for (int i = 0; i < d->nind; i++)
    d->modLabel(i) = json["modLabel"]["real"].array_items()[i].int_value();
  printf("# modLabel\n");

  d->delt_res = VectorXcd::Zero(d->nmodmax);

  d->coef_abs = VectorXcd::Zero(d->nind);
  d->expn = VectorXcd::Zero(d->nind);
  for (int i = 0; i < d->nmodmax; i++)
    d->delt_res(i) = json["delt_res"]["real"].array_items()[i].number_value() + complex<double>(1i) * json["delt_res"]["imag"].array_items()[i].number_value();
  printf("# delt_res\n");

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
  printf("# coef_abs\n");

  d->qmd1 = nNNmat(d->nmodmax);
  for (int i = 0; i < d->nmodmax; i++)
    for (int j = 0; j < NSYS; j++)
      for (int k = 0; k < NSYS; k++)
        d->qmd1[i](j, k) = json["qmd1"]["real"].array_items()[INDEX3(i, j, k)].number_value() + complex<double>(1i) * json["qmd1"]["imag"].array_items()[INDEX3(i, j, k)].number_value();
  // ##  init list/array things done ##

  // ##  init keys and gams ##
  XXimat_r keys(d->nmax, d->nind);
  XXimat_r keys1(d->nmax, d->nind);
  for (int iado = 0; iado < d->nmax; iado++) {
    for (int mp = 0; mp < d->nind; mp++) {
      keys(iado, mp) = -1;
      keys1(iado, mp) = -1;
    }
  }

  keys(0, 0) = 0;

  VectorXi key_tmp(d->nind);
  VectorXi key_tmp1(d->nind);
  for (int key_index = 0; key_index < d->nind; key_index++)
    key_tmp(key_index) = -1;

  XXimat_r g_index_list_pos(d->nmax, (2 * d->nind + 1));
  for (int iado = 0; iado < d->nmax; iado++)
    for (int mp = 0; mp < 2 * d->nind + 1; mp++)
      g_index_list_pos(iado, mp) = 0;
  // ##  init keys and gams done  ##

  printf("# keys and gams\n");

  // ##  init qmdt  ##
  d->hamt = MatrixXcd::Zero(NSYS, NSYS);
  d->qmdt = nNNmat(d->nind);

  d->qmdr = nNNmat(d->nind);
  d->qmdl = nNNmat(d->nind);
  // ##  init qmdt done ##

  // ##  init ddos  ##
  nNNmat ddos(d->nmax);
  nNNmat ddos1(d->nmax);
  nNNmat ddos2(d->nmax);
  nNNmat ddos3(d->nmax);

  printf("# ddos\n");

  vector<int> nsave(7);

  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    ddos[iado_].setZero();
  }
  ddos[0](0, 0) = 1;

  // ##  init ddos done ##

  MatrixXcd tmp, deom_mat_1;
  tmp = MatrixXcd::Zero(NSYS, NSYS);
  deom_mat_1 = MatrixXcd::Identity(NSYS, NSYS);
  // ## control init##
  const double dt = json["dt"].number_value();
  double ti = json["ti"].number_value();
  const double tf = json["tf"].number_value();
  const int nt_i = ceil(abs(ti) / dt);
  const int nt_f = ceil(abs(tf) / dt);
  const int nt = nt_i + nt_f;
  ti = -nt_i * dt;
  const int n_threads = omp_get_max_threads();
  int i;
  int iado_;
  const double dt2 = dt / 2.0;
  const double dt3 = dt / 3.0;
  const double dt4 = dt / 4.0;
  const double dt6 = dt / 6.0;
  printf("dt:%e\n", dt);
  printf("nt:%d\n", nt);
  printf("n_threads:%d\n", n_threads);
  // ## control init done##
  printf("===   INIT DONE   ===\n");
  Trie *tree = new Trie(10 * d->nind);
  // Trie tree(d->comb_list[d->combmax * (d->lmax + d->nind) + d->lmax]);
  tree->try_insert(1, 0);
  d->nddo = 1;

  // stochastic
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  cout << seed;
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0.0, 1 / sqrt(dt));
  const double alp1 = json["alp1"].number_value();
  const double alp2 = json["alp2"].number_value();
  const bool girsanov = json["girsanov"].bool_value();
  printf("alp1:%f\n", alp1);
  printf("alp2:%f\n", alp2);
  const char *frho_name = json["frho"].string_value().c_str();
  complex<double> trace(0);
  FILE *frho = fopen(frho_name, "w");
  FILE *fpolar = fopen("prop-pol.dat", "w");
  complex<double> trace_dot;
  complex<double> wx1t = 0;
  complex<double> wy1t = 0;
  complex<double> wx2t = 0;
  complex<double> wy2t = 0;
// ##final init##
#pragma omp parallel default(shared)
  for (int ii = 0; ii < nt; ii++) {
#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < d->nddo; iado_++) {
      for (int mp = 0; mp < 2 * d->nind + 1; mp++)
        g_index_list_pos(iado_, mp) = 0;
      for (int mp = 0; mp < d->nind; mp++)
        keys1(iado_, mp) = -1;
      ddos1[iado_].setZero();
    }

#pragma omp single
    {
      nsave[6] = d->nddo;
      tree->clear();
      delete tree;
      Trie *tree = new Trie(min((ullint)(d->nind * d->nind * d->nddo), d->nmax));
      d->nddo = 0;
      nsave[0] = 0;
      double x1 = distribution(generator);
      double y1 = distribution(generator);
      double x2 = distribution(generator);
      double y2 = distribution(generator);

      complex<double> z1;
      complex<double> z2;
      if (girsanov) {
        complex<double> bx1 = (double)sgn(x1) * sqrt(x1 * x1 + wx1t * wx1t) - wx1t;
        complex<double> by1 = (double)sgn(y1) * sqrt(y1 * y1 + wy1t * wy1t) - wy1t;
        complex<double> bx2 = (double)sgn(x2) * sqrt(x2 * x2 + wx2t * wx2t) - wx2t;
        complex<double> by2 = (double)sgn(y2) * sqrt(y2 * y2 + wy2t * wy2t) - wy2t;

        z1 = bx1 + (complex<double>)1i * by1;
        z2 = bx2 + (complex<double>)1i * by2;
      } else {
        z1 = x1 + (complex<double>)1i * y1;
        z2 = x2 + (complex<double>)1i * y2;
      }

      d->hamt = d->ham1;

      for (int mp = 0; mp < d->nind; mp++) {
        int m1 = d->modLabel(mp);
        d->qmdl[mp] = (alp1 + z1 * sqrt((complex<double>)alp2 / 2.)) * d->qmd1[m1] * d->qmd1[m1] + (complex<double>)1i * conj(z1) * sqrt((complex<double>)alp2 / 2.) * d->qmd1[m1];
        d->qmdr[mp] = (alp1 - z2 * sqrt((complex<double>)alp2 / 2.)) * d->qmd1[m1] * d->qmd1[m1] + (complex<double>)1i * conj(z2) * sqrt((complex<double>)alp2 / 2.) * d->qmd1[m1];
      }
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[6]; i++)
      filter_p(d, ddos, keys, ddos1, keys1, tree, i);

#pragma omp for private(iado_) schedule(dynamic, 64)
    for (iado_ = 0; iado_ < nsave[6]; iado_++) {
      keys.row(iado_) = keys1.row(iado_);
      ddos[iado_].noalias() = ddos1[iado_];
    }

#pragma omp single
    {
      nsave[1] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[1]; i++)
      construct_Mapping_p(d, keys, g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[2] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[2]; i++)
      rem_cal(ddos1, ddos, d, i, keys, g_index_list_pos);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++)
      ddos3[i].noalias() = ddos[i] + ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[2]; i++)
      ddos3[i].noalias() = ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[2]; i++)
      construct_Mapping_p(ddos3, d, keys, g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[3] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[3]; i++)
      rem_cal(ddos2, ddos3, d, i, keys, g_index_list_pos);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      ddos1[i].noalias() += ddos2[i] * 2.0;
      ddos3[i].noalias() = ddos[i] + ddos2[i] * dt2;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[2]; i++) {
      ddos1[i].noalias() += ddos2[i] * 2.0;
      ddos3[i].noalias() = ddos2[i] * dt2;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[2]; i < nsave[3]; i++) {
      ddos1[i].noalias() = ddos2[i] * 2.0;
      ddos3[i].noalias() = ddos2[i] * dt2;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[3]; i++)
      construct_Mapping_p(ddos3, d, keys, g_index_list_pos, tree, i);

#pragma omp single
    {
      nsave[4] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[4]; i++)
      rem_cal(ddos2, ddos3, d, i, keys, g_index_list_pos);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      ddos1[i].noalias() += ddos2[i] * 2.0;
      ddos3[i].noalias() = ddos[i] + ddos2[i] * dt;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[3]; i++) {
      ddos1[i].noalias() += ddos2[i] * 2.0;
      ddos3[i].noalias() = ddos2[i] * dt;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[3]; i < nsave[4]; i++)
      for (int j = 0; j < NSYS; j++)
        for (int k = 0; k < NSYS; k++) {
          ddos1[i].noalias() = ddos2[i] * 2.0;
          ddos3[i].noalias() = ddos2[i] * dt;
        }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nsave[0]; i < nsave[4]; i++)
      construct_Mapping_p(ddos3, d, keys, g_index_list_pos, tree, i);

#pragma omp single
    {

      nsave[5] = d->nddo;
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < nsave[5]; i++)
      rem_cal(ddos2, ddos3, d, i, keys, g_index_list_pos);

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave[1]; i++) {
      ddos1[i].noalias() += ddos2[i];
      ddos[i].noalias() += ddos1[i] * dt6;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[1]; i < nsave[4]; i++) {
      ddos1[i].noalias() += ddos2[i];
      ddos[i].noalias() = ddos1[i] * dt6;
    }

#pragma omp for private(i) schedule(dynamic, 64)
    for (i = nsave[4]; i < nsave[5]; i++) {
      ddos1[i].noalias() = ddos2[i];
      ddos[i].noalias() = ddos1[i] * dt6;
    }

#pragma omp single
    {
      int pos;
      if (girsanov) {
        tmp.setZero();
        for (int mp = 0; mp < d->nind; mp++) {
          gen_key(key_tmp, key_tmp1, mp, 1, d);
          int m1 = d->modLabel(mp);
          pos = tree->find(hash_bad(key_tmp1, d));
          if (pos > -1) {
            tmp += d->qmd1[m1] * ddos[pos];
          }
        }

        int pos = tree->find(1);
        trace = ddos[pos].trace();

        trace_dot = tmp.trace();
        wx1t = trace_dot / trace * sqrt((complex<double>)alp2 / 2.) * (1. - (complex<double>)1i);
        wy1t = trace_dot / trace * sqrt((complex<double>)alp2 / 2.) * (1. - (complex<double>)1i);
        wx2t = trace_dot / trace * sqrt((complex<double>)alp2 / 2.) * (1. + (complex<double>)1i);
        wy2t = -trace_dot / trace * sqrt((complex<double>)alp2 / 2.) * (1. + (complex<double>)1i);
      }
      // fprintf(fpolar, "%16.6e", ti + ii * dt);
      // fprintf(fpolar, "%20.10e%20.10e", real(trace), imag(trace));
      // fprintf(fpolar, "\n");

      // printf("%16.6e\t", ti + ii * dt);

      if (ii % (nt / 1000) == 0) {
        printf("%2.0f %% done\n", float(ii / (nt / 1000)));
        printf("nddo:%d\t", nsave[1]);
        printf("nddo:%d\t", nsave[2]);
        printf("nddo:%d\t", nsave[3]);
        printf("nddo:%d\t", nsave[4]);
        printf("nddo:%d\n", nsave[5]);
      }
      pos = tree->find(1);

      fprintf(frho, "%16.6e\t", ti + (ii + 1) * dt);
      for (int i = 0; i < NSYS; i++) {
        for (int j = 0; j < NSYS; j++) {
          // printf("%16.6e\t", ddos[pos](i, i).real());
          fprintf(frho, "%14.12e\t", ddos[pos](i, j).real());
          fprintf(frho, "%14.12e\t", ddos[pos](i, j).imag());
        }
      }
      fprintf(frho, "\n");
    }
  }
}
//

int main(int argc, char *argv[]) {
  ProfilerStart("test.prof");
  ifstream jsonFile("input.json");
  stringstream strStream;
  strStream << jsonFile.rdbuf();
  string jsonStr = strStream.str();
  string err;
  const Json json = Json::parse(jsonStr, err);
  if (!err.empty()) {
    printf("Error in parsing input file: %s\n", err.c_str());
    return 1;
  }

  init_deom(json, argc, argv);
  ProfilerStop();
}