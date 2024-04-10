/* assert */
#include <assert.h>
/* errno */
#include <errno.h>
/* fopen, fscanf, fprintf, fclose */
#include <stdio.h>
/* EXIT_SUCCESS, EXIT_FAILURE, malloc, free */
#include <stdlib.h>
/*imports omp*/
#include <omp.h>
/*imports time*/
#include <time.h>


static int create_mat(size_t const nrows, size_t const ncols, double** const matp) {
    double* mat = NULL;
    if (!(mat = (double*)malloc(nrows * ncols * sizeof(*mat)))) {
        goto cleanup;
    }
    /** Initialize matrix with random values **/
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            mat[(i * ncols) + j] = (double)(rand() % 1000) / 353.0;
        }
    }
    /** End random initialization **/
    *matp = mat;
    return 0;
    cleanup:
    free(mat);
    return -1;
}

//nrows = n, ncols = m, ncols2 = t, mat1 = A, mat 2 = B, return mat = C
static int mult_mat_tile(size_t const n, size_t const m, size_t const p, double const* const A, double const* const B, double** const Cp, int numThreads) {
    //define variables
    size_t i, j, k, t1, t2;
    int maxThreads = omp_get_max_threads();
    double start; 
    double end; 
    double sum;
    int tileSize;
    double* C = NULL;
    omp_set_num_threads(numThreads);


    //allocate space required
    if (!(C = (double*)malloc(n * p * sizeof(*C)))) {
        goto cleanup;
    }
    
    tileSize = 45;
    
    //parallel for loop to multiply tiled matrices
    #pragma omp parallel for shared(C, sum) private(t1, t2, i, j, k)
    for(t1 = 0; t1 < n; t1 += tileSize) {
        for(t2 = 0; t2 < p; t2 += tileSize) {
            for (i = t1; i < n && i < t1 + tileSize; ++i) {
                for(int j = t2; j < p && j < t2 + tileSize; ++j) {
                    for (k = 0, sum = 0.0; k < m; ++k) {
                        sum += A[i * m + k] * B[k * p + j];
                    }
                    C[i * p + j] = sum;
                }
            }
        }
    }


    
    *Cp = C;
    return 0;
cleanup:
    free(C);
    /*failure:*/
    return -1;
}


//nrows = n, ncols = m, ncols2 = t, mat1 = A, mat 2 = B, return mat = C
static int mult_mat(size_t const n, size_t const m, size_t const p, double const* const A, double const* const B, double** const Cp, int numThreads) {
    //define variables
    size_t i, j, k;
    int maxThreads = omp_get_max_threads();
    double start; 
    double end; 
    double sum;
    double* C = NULL;
    omp_set_num_threads(numThreads);


    //allocate space required
    if (!(C = (double*)malloc(n * p * sizeof(*C)))) {
        goto cleanup;
    }
    
    //parallel for loop to multiply matrices
    #pragma omp parallel for shared(C, sum) private(i, j, k)
    for (i = 0; i < n; ++i) {
        for (j = 0; j < p; ++j) {
            for (k = 0, sum = 0.0; k < m; ++k) {
                sum += A[i * m + k] * B[k * p + j];
            }
            C[i * p + j] = sum;
        }
    }

    

    *Cp = C;
    return 0;
cleanup:
    free(C);
    /*failure:*/
    return -1;
}




int main(int argc, char* argv[]) {
    // size_t stored an unsigned integer
    size_t nrows, ncols, ncols2;
    double *A = NULL, *B = NULL, *C = NULL;
    int k = 1;
    if (argc != 4) {
        fprintf(stderr, "usage: matmult nrows ncols ncols2\n");
        goto failure;
    }
    nrows = atoi(argv[1]);
    ncols = atoi(argv[2]);
    ncols2 = atoi(argv[3]);

    double start; 
    double end; 


    //2D parallelization
    if (create_mat(nrows, ncols, &A)) {
        perror("error");
        goto failure;
    }
    if (create_mat(ncols, ncols2, &B)) {
        perror("error");
        goto failure;
    }

    int numThreads = 16;
    
    
    for(int i = 1; i < 29; i++) {
        if(i % 4 == 0 || i == 1 || i == 2) {
            numThreads = i;
            //tiles
            start = omp_get_wtime(); 
            if(mult_mat_tile(nrows, ncols, ncols2, A, B, &C, numThreads)) {
                perror("error");
                goto failure;
            }
            end = omp_get_wtime(); 
            printf("Total time tiled = %f\n", end - start);

            //normal
            start = omp_get_wtime();
            if(mult_mat(nrows, ncols, ncols2, A, B, &C, numThreads)) {
                perror("error");
                goto failure;
            }
            end = omp_get_wtime(); 
            printf("Total time not tiled= %f\n", end - start);
            printf("Number of Threads = %d\n", numThreads);
        }
    }
    


    free(A);
    free(B);
    free(C);
    return EXIT_SUCCESS;
failure:
    if (A) {
        free(A);
    }
    if (B) {
        free(B);
    }
    if (C) {
        free(C);
    }
    return EXIT_FAILURE;
}
