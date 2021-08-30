#pragma once

#include "deom.hpp"
#include "algebra.hpp"
#include <Eigen/Eigenvalues>
#include <experimental/random>

// syl： solve syl equation AX-XB=C
// note here A and B are both hermmitian, so it canbe deocomposite to URU^T VSV^T where R and S are diag matrices, the equation then canbe convert to URU^TX-XVSV^T=C multi U^T on left and V on right so we can get RU^TXV-U^TXVS=U^TCV, denote U^TXV as Y and U^TCV as D,so \sum_j (R_ij Y_jk - Y_ij S_jk) = D_ik and note that R and S are diag matrices so R_ii Y_ik - Y_ik S_kk = D_ik.
// SelfAdjointEigenSolver(A,B) -> R S U V
// Y_ik = (U^TCV)_ik / (R_ii - S_kk)
// X = UYV^T
// check AX-XB=C
// note that syl_complex will deal with no SelfAdjoint situation, but we will no use that in normal scenario.
// pref
// syl :1Gflops 10000000 at 3.6s
// syl_debug :1Gflops 10000000 at 3.85s
// syl_complex :2Gflops 10000000 at 5.4s

void syl(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  SelfAdjointEigenSolver<MatrixNcd> eigensolverA(A);
  SelfAdjointEigenSolver<MatrixNcd> eigensolverB(B);
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = eigensolverA.eigenvectors();
  MatrixNcd V = eigensolverB.eigenvectors();
  MatrixNcd Y = MatrixNcd::Zero();
  MatrixNcd D = U.transpose() * C * V;

  for (int i = 0; i < NSYS; i++) {
    for (int j = 0; j < NSYS; j++) {
      Y(i, j) = D(i, j) / (R(i) - S(j));
    }
  }
  X = U * Y * V.transpose();
}

void syl_debug(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  SelfAdjointEigenSolver<MatrixNcd> eigensolverA(A);
  SelfAdjointEigenSolver<MatrixNcd> eigensolverB(B);
  if (eigensolverA.info() != Success || eigensolverB.info() != Success)
    abort();
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = eigensolverA.eigenvectors();
  MatrixNcd V = eigensolverB.eigenvectors();
  MatrixNcd Y = MatrixNcd::Zero();
  MatrixNcd D = U.transpose() * C * V;

  for (int i = 0; i < NSYS; i++) {
    for (int j = 0; j < NSYS; j++) {
      Y(i, j) = D(i, j) / (R(i) - S(j));
    }
  }
  X = U * Y * V.transpose();
  if (is_valid(A * X - X * B - C)) {
    exit(1);
  }
}

void syl_complex(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  ComplexEigenSolver<MatrixNcd> eigensolverA(A);
  ComplexEigenSolver<MatrixNcd> eigensolverB(B);
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = eigensolverA.eigenvectors();
  MatrixNcd V = eigensolverB.eigenvectors();
  MatrixNcd Y = MatrixNcd::Zero();
  MatrixNcd D = U.inverse() * C * V;

  for (int i = 0; i < NSYS; i++) {
    for (int j = 0; j < NSYS; j++) {
      Y(i, j) = D(i, j) / (R(i) - S(j));
    }
  }
  X = U * Y * V.inverse();
}

void syl_complex_debug(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  ComplexEigenSolver<MatrixNcd> eigensolverA(A);
  ComplexEigenSolver<MatrixNcd> eigensolverB(B);
  if (eigensolverA.info() != Success || eigensolverB.info() != Success)
    abort();
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = eigensolverA.eigenvectors();
  MatrixNcd V = eigensolverB.eigenvectors();
  MatrixNcd Y = MatrixNcd::Zero();
  MatrixNcd D = U.inverse() * C * V;

  for (int i = 0; i < NSYS; i++) {
    for (int j = 0; j < NSYS; j++) {
      Y(i, j) = D(i, j) / (R(i) - S(j));
    }
  }
  X = U * Y * V.inverse();
  if (is_valid(A * X - X * B - C)) {
    exit(1);
  }
}

// int main(int argc, const char **argv) {
//   MatrixNcd A;
//   MatrixNcd B;
//   MatrixNcd C;
//   MatrixNcd X;
//   std::srand(std::time(0));
//   A(0, 0) = std::experimental::randint(-999, 999);
//   A(1, 1) = std::experimental::randint(-999, 999);
//   A(0, 1) = std::experimental::randint(-999, 999);
//   A(1, 0) = A(0, 1);

//   B(0, 0) = std::experimental::randint(-999, 999);
//   B(1, 1) = std::experimental::randint(-999, 999);
//   B(0, 1) = std::experimental::randint(-999, 999);
//   B(1, 0) = B(0, 1);

//   C(0, 0) = (double)std::experimental::randint(-999, 999) + (complex<double>)1i * (double)std::experimental::randint(-999, 999);
//   C(1, 1) = (double)std::experimental::randint(-999, 999) + (complex<double>)1i * (double)std::experimental::randint(-999, 999);
//   C(0, 1) = (double)std::experimental::randint(-999, 999) + (complex<double>)1i * (double)std::experimental::randint(-999, 999);
//   C(1, 0) = (double)std::experimental::randint(-999, 999) + (complex<double>)1i * (double)std::experimental::randint(-999, 999);
//   for (size_t i = 0; i < 10000000; i++) {
//     syl_complex(X, A, B, C);
//   }
//   return 0;
// }
