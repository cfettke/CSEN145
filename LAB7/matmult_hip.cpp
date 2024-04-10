#include <iostream>
#include <hip/hip_runtime.h>
#include <stdexcept>
#include "matmult.h"

#define BLOCK_SIZE 16

inline void hipAssert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess) 
   {
      fprintf(stderr,"HIP GPU assert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) {
         exit(code);
      }
   }
}

#define hipErrorCheck(ans) { hipAssert((ans), __FILE__, __LINE__); }

/**
 * Matrix multiplication kernel
 * @param a GPU device pointer to an m X n matrix (A)
 * @param b GPU device pointer to a n X k matrix (B)
 * @param c GPU device output pointer to an m X k matrix (C) 
 * @param m nrows of A and C
 * @param n ncols of A and nrows of B
 * @param k ncols of B and C
 * 
 * Note:
 *     grid and block should be configured as:
 *         dim3 dimGrid((k + BLOCK_SIZE - 1) / BLOCK_SIZE, (m + BLOCK_SIZE - 1) / BLOCK_SIZE);
 *         dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
 */

__global__ void hip_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k)
{ 
  int row = hipBlockIdx_y * hipBlockDim_y + hipThreadIdx_y; 
  int col = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if( col < k && row < m) 
  {
    double v = 0.0;
    for(int i = 0; i < n; i++) 
    {
      v += a[row * n + i] * b[i * k + col];
    }
    c[row * k + col] = v;
  }
} 

int main(int argc, char const *argv[])
{
    
  size_t nrows, ncols, ncols2;
  bool verify = argc > 4;
  
  if (argc < 4) {
    fprintf(stderr, "usage: matmult nrows ncols ncols2\n");
    return -1;
  }

  nrows = atoi(argv[1]);
  ncols = atoi(argv[2]);
  ncols2 = atoi(argv[3]);

  std::cout << "Matrix multiplication A(" << nrows << ", " << ncols << ") x B(" << ncols << ", " << ncols2 << ")" << std::endl;
  
  size_t memsize_A = sizeof(double) * nrows * ncols;
  size_t memsize_B = sizeof(double) * ncols * ncols2;
  size_t memsize_C = sizeof(double) * nrows * ncols2;

  // create host matrices
  auto A = create_mat(nrows, ncols, true);
  auto B = create_mat(ncols, ncols2, true);
  auto C = create_mat(nrows, ncols2, false);

  // create device matrices
  double *dA, *dB, *dC;
  hipErrorCheck(hipMalloc((void **) &dA, memsize_A))
  hipErrorCheck(hipMalloc((void **) &dB, memsize_B))
  hipErrorCheck(hipMalloc((void **) &dC, memsize_C))

  // copy host memory to device
  hipErrorCheck(hipMemcpy(dA, A, memsize_A, hipMemcpyHostToDevice))
  hipErrorCheck(hipMemcpy(dB, B, memsize_B, hipMemcpyHostToDevice))

  // set up kernel execution parameters
  unsigned int grid_rows = (nrows + BLOCK_SIZE - 1) / BLOCK_SIZE;
  unsigned int grid_cols = (ncols2 + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 dimGrid(grid_cols, grid_rows);
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

  // execute kernel
  hip_matmult <<< dimGrid, dimBlock >>> (dA, dB, dC, nrows, ncols, ncols2);

  // transfer results from device to host
  hipErrorCheck(hipMemcpy(C, dC, memsize_C, hipMemcpyDeviceToHost))

  // optionally verify
  if(verify){
    // execute serial program on CPU
    auto C2 = create_mat(nrows, ncols2, false);
    serial_matmult(A, B, C2, nrows, ncols, ncols2);
    verify_result(C, C2, nrows, ncols2);
    free(C2);
  }

  // free memory
  hipErrorCheck(hipFree(dA))
  hipErrorCheck(hipFree(dB))
  hipErrorCheck(hipFree(dC))
  free(A);
  free(B);
  free(C);

  return 0;
}
