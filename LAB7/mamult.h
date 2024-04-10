#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>

//#include <openacc.h>

using namespace std;

/**
 * Matrix multiplication serial algorithm
 * 
 * @param a host device pointer to an m X n matrix (A)
 * @param b host device pointer to a n X k matrix (B)
 * @param c host device output pointer to an m X k matrix (C) 
 * @param m nrows of A and C
 * @param n ncols of A and nrows of B
 * @param k ncols of B and C
*/
void serial_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k) {
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < k; ++j)
    {
      double v = 0.0;
      for (int h = 0; h < n; ++h)
      {
        v += a[i * n + h] * b[h * k + j];
      }
      c[i * k + j] = v;
    }
  }
}


void parallel_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k) {
  #pragma omp parallel for shared(c)
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < k; ++j)
    {
      double v = 0.0;
      for (int h = 0; h < n; ++h)
      {
        v += a[i * n + h] * b[h * k + j];
      }
      c[i * k + j] = v;
    }
  }
}

    
void parallel_ACC_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k) {
  //#pragma acc parallel loop
  //#pragma acc parallel /*kernels loop*/ copyin(a[0: m*n], b[0: n*k]) copyout(a[0: m*k])
  #pragma acc data copyin(a[0: m*n], b[0: n*k]) copyout(a[0: m*k])
  {
    #pragma acc region
    {
      #pragma acc loop
      {
        for (int i = 0; i < m; ++i)
        {
          for (int j = 0; j < k; ++j)
          {
            double v = 0.0;
            for (int h = 0; h < n; ++h)
            {
              v += a[i * n + h] * b[h * k + j];
            }
            c[i * k + j] = v;
          }
        }
      }
    }
  }
}


/** 
 * Create a matrix of size m x n and optionally initialize it with values
 */
double * create_mat(unsigned int m, unsigned int n, bool init=false)
{
  double * mat = NULL;
  if (!(mat = (double*) malloc(m*n*sizeof(*mat)))){
    throw std::runtime_error("Could not allocate matrix.");
  }
  if(init){
    for(size_t i = 0; i < m; i++){
      for (size_t j = 0; j < n; j++){
        mat[(i * n) + j] = (double)(i*n+j);
      }
    }
  }
  return mat;
}

/**
 * Verify if matrices C and C2 contain the same values
*/
void verify_result(double *C, double *C2, unsigned int m, unsigned int n){
  for(size_t i = 0; i < m; i++){
    for (size_t j = 0; j < n; j++){
      if( abs(C[i * n + j] - C2[i * n + j]) > 1e-5 ){
        cout << "Matrices C and C2 do not match!!!" << endl;
        return;
      }
    }
  }
  cout << "Matrices C and C2 match." << endl;
}