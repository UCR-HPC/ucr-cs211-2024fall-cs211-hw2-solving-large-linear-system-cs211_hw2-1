#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h> // OpenMP for parallelization

#define BLOCK_SIZE 128

void mydgemm(int n, double *A, double *B, double *C) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i += BLOCK_SIZE) {
        for (int j = 0; j < n; j += BLOCK_SIZE) {
            for (int k = 0; k < n; k += BLOCK_SIZE) {
                for (int ii = i; ii < i + BLOCK_SIZE && ii < n; ++ii) {
                    for (int jj = j; jj < j + BLOCK_SIZE && jj < n; ++jj) {
                        double sum = 0.0;
                        for (int kk = k; kk < k + BLOCK_SIZE && kk < n; ++kk) {
                            sum += A[ii * n + kk] * B[kk * n + jj];
                        }
                        C[ii * n + jj] += sum;
                    }
                }
            }
        }
    }
}

void mydgetrf_block(int n, double *A) {
    for (int k = 0; k < n; k += BLOCK_SIZE) {
        int end = k + BLOCK_SIZE < n ? k + BLOCK_SIZE : n;

        // Factorize the diagonal block
        for (int i = k; i < end; ++i) {
            for (int j = k; j < i; ++j) {
                A[i * n + j] /= A[j * n + j];
                for (int l = j + 1; l < end; ++l) {
                    A[i * n + l] -= A[i * n + j] * A[j * n + l];
                }
            }
        }

        // Update the trailing submatrix
        #pragma omp parallel for collapse(2)
        for (int i = end; i < n; ++i) {
            for (int j = k; j < end; ++j) {
                A[i * n + j] /= A[j * n + j];
                for (int l = end; l < n; ++l) {
                    A[i * n + l] -= A[i * n + j] * A[j * n + l];
                }
            }
        }

        // Perform matrix multiplication for the trailing submatrix
        #pragma omp parallel for collapse(2)
        for (int i = end; i < n; i += BLOCK_SIZE) {
            for (int j = end; j < n; j += BLOCK_SIZE) {
                for (int ii = i; ii < i + BLOCK_SIZE && ii < n; ++ii) {
                    for (int jj = j; jj < j + BLOCK_SIZE && jj < n; ++jj) {
                        double sum = 0.0;
                        for (int kk = k; kk < end; ++kk) {
                            sum += A[ii * n + kk] * A[kk * n + jj];
                        }
                        A[ii * n + jj] -= sum;
                    }
                }
            }
        }
    }
}

int main() {
    int n = 5000; // Example size
    FILE *pad_file = fopen("pad.txt", "r");
    int pad;
    fscanf(pad_file, "%d", &pad);
    fclose(pad_file);
    n = ((n + pad - 1) / pad) * pad;
    printf("n=%d, pad=%d\n", n, pad);

    double *A = (double *)malloc(n * n * sizeof(double));

    // Initialize matrix A with some values
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = ((double)rand() / RAND_MAX) * 2 - 1;
        }
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);

    mydgetrf_block(n, A);

    gettimeofday(&end, NULL);
    double elapsed_time = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
    printf("Time taken: %lf seconds\n", elapsed_time);

    free(A);
    return 0;
}