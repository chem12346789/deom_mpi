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
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread> // std::thread
#include <unistd.h>

#define NSYS 2

using namespace json11;
using namespace std;
using namespace Eigen;

typedef unsigned long long int ullint;
typedef long long int llint;
typedef Eigen::Matrix<complex<double>, NSYS, NSYS> MatrixNcd;
typedef Eigen::Matrix<ullint, Eigen::Dynamic, Eigen::Dynamic> nullmat;
typedef std::vector<MatrixNcd, Eigen::aligned_allocator<MatrixNcd>> nNNmat;
typedef std::vector<MatrixXcd, Eigen::aligned_allocator<MatrixXcd>> nXXmat;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXimat_r;
typedef Eigen::Matrix<uint16_t, Eigen::Dynamic, 1, Eigen::RowMajor> vetui16;
typedef std::vector<vetui16, Eigen::aligned_allocator<MatrixXcd>> XXui16mat_r;


