#include <stdio.h>
#include <stdlib.h>
#include <math.h> // 确保这个在最上面
#include <time.h>  // 确保时间函数的声明

#ifndef __MY_C__
#define __MY_C__


int mydgetrf(double *A, int *ipiv, int n) {
    for (int i = 0; i < n - 1; i++) {
        // Pivoting
        int maxind = i;
        double max = fabs(A[i * n + i]);
        printf("Step %d: Current max = %f at row %d\n", i, max, i);

        for (int t = i + 1; t < n; t++) {
            if (fabs(A[t * n + i]) > max) {
                maxind = t;
                max = fabs(A[t * n + i]);
            }
        }

        if (max == 0) {
            printf("Matrix is singular at step %d\n", i);
            return 0; // 矩阵是不可约的
        } else {
            if (maxind != i) {
                // 交换行
                printf("Swapping rows %d and %d\n", i, maxind);
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

        // LU 分解
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];
            printf("Row %d after division: ", j);
            for (int k = 0; k < n; k++) {
                printf("%f ", A[j * n + k]);
            }
            printf("\n");

            for (int k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }

        printf("Matrix A after step %d:\n", i);
        for (int m = 0; m < n; m++) {
            for (int k = 0; k < n; k++) {
                printf("%f ", A[m * n + k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    return 1; // 结果有效
}


void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) {
    if (UPLO == 'L') {
        double *tempB = (double *)malloc(n * sizeof(double)); // 创建一个新数组
        for (int i = 0; i < n; i++) {
            tempB[i] = B[ipiv[i]]; // 根据 ipiv 重排
            printf("tempB[%d] after pivoting: %f\n", i, tempB[i]);
        }
        for (int i = 0; i < n; i++) {
            B[i]=tempB[i];
            printf("B[%d] before forward substitution (after pivoting): %f\n", i, B[i]);
        }
        
        // 前向替换
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                B[i] -= A[i * n + j] * B[j];
                printf("B[i]: %f\n",  B[i]);
                printf("A[i * n + j]]: %f\n", A[i * n + j]);
            }
            printf("B[%d] after forward substitution: %f\n", i, B[i]);
        }
    }else if (UPLO == 'U') {
        // 后向替换
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i * n + j] * B[j];
            }
            if (fabs(A[i * n + i]) < 1e-10) {
                printf("Error: Division by zero at row %d\n", i);
                return;
            }
            printf("B[i]: %f\n",  B[i]);
            printf("A[i * n + i]]: %f\n", A[i * n + i]);
            B[i] = (B[i] - sum) / A[i * n + i];
            printf("#########B[%d] after backward substitution: %f\n", i, B[i]);
        }
        // double *tempB = (double *)malloc(n * sizeof(double)); // 创建一个新数组
        // for (int i = 0; i < n; i++) {
        //     tempB[ipiv[i]] = B[i];
        //     printf("tempB[%d] after pivoting: %f\n", i, tempB[i]);
        // }
        // for (int i = 0; i < n; i++) {
        //     B[i]=tempB[i];
        //     printf("B[%d] after backword substitution (after pivoting): %f\n", i, B[i]);
        // }
    }
}


void my_f(double *A, double *B, int n) {
    int *ipiv = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        ipiv[i] = i;
    }

    printf("Initial Matrix A:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i * n + j]);
        }
        printf("\n");
    }

    printf("Initial Vector B:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }
    printf("\n\n");

    if (mydgetrf(A, ipiv, n) == 0) {
        printf("LU factorization failed: coefficient matrix is singular.\n");
        free(ipiv);
        return;
    }

    printf("Matrix A after LU factorization:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i * n + j]);
        }
        printf("\n");
    }
    
    printf("Pivot indices:\n");
    for (int i = 0; i < n; i++) {
        printf("%d ", ipiv[i]);
    }
    printf("\n\n");

    mydtrsv('L', A, B, n, ipiv); // 前向替换
    printf("Vector B after forward substitution:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }
    printf("\n\n");

    mydtrsv('U', A, B, n, ipiv); // 后向替换
    printf("Vector B after backward substitution (Solution X):\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }
    printf("\n\n");

    free(ipiv);
}

#endif


int main() {
    int n = 3; 
    double *A = (double *)malloc(n * n * sizeof(double));
    double *B = (double *)malloc(n * sizeof(double));
    

    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = (rand() % 10) + 1; 
        }
        B[i] = (rand() % 10) + 1; 
    }

    printf("Matrix A:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i * n + j]);
        }
        printf("\n");
    }

    printf("Vector B:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }
    printf("\n");

    my_f(A, B, n); 

    printf("Solution X:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }
    printf("\n");

    free(A);
    free(B);

    return 0;
}
