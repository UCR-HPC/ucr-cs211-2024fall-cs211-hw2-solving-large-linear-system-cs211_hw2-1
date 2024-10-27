#ifndef __MY_C__
#define __MY_C__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int mydgetrf(double *A, int *ipiv, int n) {
    for (int i = 0; i < n - 1; i++) {
        // Pivoting
        int maxind = i;
        double max = fabs(A[i * n + i]);
        //printf("Step %d: Current max = %f at row %d\n", i, max, i);

        for (int t = i + 1; t < n; t++) {
            if (fabs(A[t * n + i]) > max) {
                maxind = t;
                max = fabs(A[t * n + i]);
            }
        }

        if (max == 0) {
            printf("Matrix is singular at step %d\n", i);
            return 0; 
        } else {
            if (maxind != i) {

                //printf("Swapping rows %d and %d\n", i, maxind);
                int temp = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temp;

                for (int j = 0; j < n; j++) {
                    double temp_value = A[i * n + j];
                    A[i * n + j] = A[maxind * n + j];
                    A[maxind * n + j] = temp_value;
                }
            }
        }

        // LU 
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];
            for (int k = 0; k < n; k++) {
            }

            for (int k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    return 1;
}


void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) {
    if (UPLO == 'L') {
        double *tempB = (double *)malloc(n * sizeof(double)); 
        for (int i = 0; i < n; i++) {
            tempB[i] = B[ipiv[i]]; 
        }
        for (int i = 0; i < n; i++) {
            B[i]=tempB[i];
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                B[i] -= A[i * n + j] * B[j];
            }
        }
    }else if (UPLO == 'U') {
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i * n + j] * B[j];
            }
            if (fabs(A[i * n + i]) < 1e-10) {
                printf("Error: Division by zero at row %d\n", i);
                return;
            }
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
}


void my_f(double *A, double *B, int n) {
    int *ipiv = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        ipiv[i] = i;
    }

    if (mydgetrf(A, ipiv, n) == 0) {
        printf("LU factorization failed: coefficient matrix is singular.\n");
        free(ipiv);
        return;
    }

    mydtrsv('L', A, B, n, ipiv);
    mydtrsv('U', A, B, n, ipiv); 

    free(ipiv);
}

#endif