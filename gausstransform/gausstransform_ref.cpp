/* 
 * This file contains a few functions from:
 *
 * Robust Point Set Registration Using Mixture of Gaussians
 * Copyright (C) 2008, Bing Jian and Baba C. Vemuri
 * https://github.com/bing-jian/gmmreg
 *
 * The full algorithm is described in the following ICCV'05 paper:
 *      A Robust Algorithm for Point Set Registration Using Mixture of Gaussians,
 *      Bing Jian and Baba C. Vemuri,
 *      10th IEEE International Conference on Computer Vision (ICCV 2005),
 *      17-20 October 2005, Beijing, China, pp. 1246-1251.
 *
 * The code in this file is licensed under
 *      the GNU General Public License v3.0 (GPLv3)
 */

#define SQR(X)  ((X)*(X))

template<typename T>
T GaussTransform(const T* A, const T* B,
    int m, int n, int dim, T scale) {
  T cross_term = 0;
  scale = SQR(scale);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      T dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij += SQR(A[i * dim + d] - B[j * dim + d]);
      }
      T cost_ij = exp(-1.0 * dist_ij / scale);
      cross_term += cost_ij;
    }
    /* printf("cross_term = %.3f\n", cross_term);  */
  }
  return cross_term / (m * n);
}

template<typename T>
T GaussTransform(const T* A, const T* B,
    int m, int n, int dim, T scale, T* grad) {
  T cross_term = 0;
  for (int i = 0; i < m * dim; ++i) {
    grad[i] = 0;
  }

  scale = SQR(scale);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      T dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij +=  SQR(A[i * dim + d] - B[j * dim + d]);
      }
      T cost_ij = exp(-1.0 * dist_ij / scale);
      for (int d = 0; d < dim; ++d){
        grad[i * dim + d] -= cost_ij * 2.0 * (A[i * dim + d] - B[j * dim + d]);
      }
      cross_term += cost_ij;
    }
    /* printf("cross_term = %.3f\n", cross_term);  */
  }
  scale *= m * n;
  for (int i = 0; i < m * dim; ++i) {
    grad[i] /= scale;
  }
  return cross_term / (m * n);
}

