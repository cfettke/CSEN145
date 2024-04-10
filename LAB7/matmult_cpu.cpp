#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>
#include "matmult.h"

int main(int argc, char const *argv[])
{
    
  size_t nrows, ncols, ncols2;
  
  if (argc < 4) {
    fprintf(stderr, "usage: matmult nrows ncols ncols2\n");
    return -1;
  }

  nrows = atoi(argv[1]);
  ncols = atoi(argv[2]);
  ncols2 = atoi(argv[3]);

  cout << "Matrix multiplication A(" << nrows << ", " << ncols << ") x B(" << ncols << ", " << ncols2 << ")" << endl;

  // create host matrices
  auto A = create_mat(nrows, ncols, true);
  auto B = create_mat(ncols, ncols2, true);
  auto C = create_mat(nrows, ncols2, false);

  serial_matmult(A, B, C, nrows, ncols, ncols2);

  // free memory
  free(A);
  free(B);
  free(C);

  return 0;
}