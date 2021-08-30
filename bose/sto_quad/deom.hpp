#ifndef DEOM_H_
#define DEOM_H_

#include "hashmap.hpp"
#include "json11.hpp"
#include <Eigen/Dense>
#include <atomic>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using Eigen::Matrix2cd;
using Eigen::Matrix3cd;
using Eigen::Matrix4cd;
using Eigen::MatrixXcd;
using Eigen::MatrixXi;
using Eigen::VectorXcd;
using Eigen::VectorXi;

#define NSYS 2

using namespace json11;
using namespace std;

typedef unsigned long long int ullint;
typedef long long int llint;
typedef Eigen::Matrix<complex<double>, NSYS, NSYS> MatrixNcd;
typedef Eigen::Matrix<ullint, Eigen::Dynamic, Eigen::Dynamic> nullmat;
typedef std::vector<MatrixNcd, Eigen::aligned_allocator<MatrixNcd>> nNNmat;
typedef std::vector<MatrixXcd, Eigen::aligned_allocator<MatrixXcd>> nXXmat;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXimat_r;
typedef struct {
  ullint nmax;
  int nsys;
  int nind;
  int nddo;
  float ferr;
  int combmax;
  int lmax;
  int nmod;
  int nmodmax;
  int nham;
  nullmat comb_list;
  VectorXi modLabel;
  mutex nddo_mutex;
  Matrix2cd ham1;
  Matrix2cd hamt;
  VectorXcd delt_res;
  nNNmat qmd1;
  nNNmat qmdt;
  nNNmat qmdl;
  nNNmat qmdl_s;
  nNNmat qmdr;
  nNNmat qmdr_s;
  VectorXcd coef_abs;
  nXXmat coef_lft;
  nXXmat coef_rht;
  VectorXcd expn;
  nNNmat sdip;
  nNNmat pdip;
  VectorXcd bdip;
} DEOM;

// Please note that due to the convenience of inspection, we always keep the non-parallel version.
#endif