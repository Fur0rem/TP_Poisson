/**
 * @file test_facto_LU.c
 * @brief Test factorisation LU
 */

#include "atlas_headers.h"
#include "lib_poisson1D.h"
#include "test.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
        int nb_rows = 3;
        int la = 3;
        int kv = 1;
        int lab = nb_rows + kv;
        double* poisson = (double*)malloc(lab * la * sizeof(double));
        set_GB_operator_colMajor_poisson1D(poisson, &lab, &la, &kv);
        print_mat_col_major(poisson, lab, la);

        // LU factorisation
        int* ipiv = (int*)malloc(la * sizeof(int));
        int info;
        dgbtrftridiag(&la, &la, &kv, &kv, poisson, &lab, ipiv, &info);
        ASSERT(info == 0);
        print_mat_col_major(poisson, lab, la);

        // Extract L and U
        double* L = (double*)malloc(la * la * sizeof(double));
        double* U = (double*)malloc(la * la * sizeof(double));

        // All elements on the diagonal of L are 1
        for (int i = 0; i < la; i++) {
                L[i * la + i] = 1;
        }

        // All the elements of the diagonal of A are the diagonal of U
        for (int i = 0; i < la; i++) {
                int idx = indexABCol(i, i, &lab);
                U[i * la + i] = poisson[idx];
        }

        // All the elements of the upper triangle of A are the upper triangle of U
        for (int i = 0; i < la; i++) {
                for (int j = i + 1; j < la; j++) {
                        int idx = indexABCol(i, j, &lab);
                        U[j * la + i] = poisson[idx];
                }
        }

        // All the elements of the lower triangle of A are the lower triangle of L
        for (int i = 0; i < la; i++) {
                for (int j = 0; j < i; j++) {
                        int idx = indexABCol(i, j, &lab);
                        L[j * la + i] = poisson[idx];
                }
        }

        printf("L:\n");
        print_mat_row_major((double*)L, 3, 3);
        printf("U:\n");
        print_mat_row_major((double*)U, 3, 3);

        double __expected_L[3][3] = {
            {1, 0, 0},
            {-0.5, 1, 0},
            {0, -2.0 / 3.0, 1},
        };
        double* expected_L = flattened_row_major(3, 3, __expected_L);

        ASSERT(mat_equals(expected_L, L, 3, 3));

        double __expected_U[3][3] = {
            {2, -1, 0},
            {0, 1.5, -1},
            {0, 0, 1.333333},
        };
        double* expected_U = flattened_row_major(3, 3, __expected_U);

        ASSERT(mat_equals(expected_U, U, 3, 3));

        // Multiply L and U
        double* LU = (double*)malloc(la * la * sizeof(double));
        for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                                LU[i * la + j] += L[k * la + i] * U[j * la + k];
                        }
                }
        }

        printf("LU:\n");
        print_mat_row_major(LU, 3, 3);

        double __expected_LU[3][3] = {
            {2, -1, 0},
            {-1, 2, -1},
            {0, -1, 2},
        };
        double* expected_LU = flattened_row_major(3, 3, __expected_LU);

        ASSERT(mat_equals(expected_LU, LU, 3, 3));

        printf("Test factorisation tridiagonal passed\n");
}