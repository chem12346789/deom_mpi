#pragma once
#include "hashmap.hpp"
#include "json11.hpp"
#include "pulse.hpp"
// #include "dipole.hpp"
#include <Eigen/Dense>
#include <atomic>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gperftools/profiler.h>
#include <iomanip>
#include <iostream>
#include <openssl/md5.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread> // std::thread
#include <unistd.h>

#define NSYS 4

using namespace json11;
using namespace std;
using namespace Eigen;

typedef unsigned long long int ullint;
// max number of bytes used to store a pointer in 64 bit
typedef long long int llint;
typedef Eigen::Matrix<complex<double>, NSYS, NSYS> MatrixNcd;
typedef Eigen::Matrix<ullint, Eigen::Dynamic, Eigen::Dynamic> nullmat;
typedef std::vector<MatrixNcd, Eigen::aligned_allocator<MatrixNcd>> nNNmat;
typedef std::vector<MatrixXcd, Eigen::aligned_allocator<MatrixXcd>> nXXmat;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXimat_r;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXbmat_r;
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic, Eigen::RowMajor> XXbvet_r;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic, Eigen::RowMajor> XXivet_r;
// typedef Eigen::Matrix<uint16_t, Eigen::Dynamic, 1, Eigen::RowMajor> vetui16;
// typedef std::vector<vetui16, Eigen::aligned_allocator<MatrixXcd>> XXui16mat_r;


