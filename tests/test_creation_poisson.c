/**
 * @file test_creation_poisson.c
 * @brief Test creation of poisson matrix

 */

#include "atlas_headers.h"
#include "lib_poisson1D.h"
#include "test.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
        // Test de la matrice identité
        {
                int nb_rows = 1;
                int la = 3;
                int kv = 1;
                int lab = nb_rows + kv;
                double* mat = (double*)malloc(lab * la * sizeof(double));
                set_GB_operator_colMajor_poisson1D_Id(mat, &lab, &la, &kv);
                print_mat_col_major(mat, lab, la);

                double __expected[2][3] = {
                    {0, 0, 0},
                    {1, 1, 1},
                };
                double* expected = flattened_row_major(2, 3, __expected);
                ASSERT(mat_equals(mat, expected, 2, 3));
                printf("Test Creation poisson identity matrix passed\n");

                free(expected);
                free(mat);
        }

        // Test de la matrice de poisson
        {
                int nb_rows = 3;
                int la = 7;
                int kv = 2;
                int lab = nb_rows + kv;

                double* mat = (double*)malloc(lab * la * sizeof(double));
                set_GB_operator_colMajor_poisson1D(mat, &lab, &la, &kv);
                print_mat_col_major(mat, lab, la);

                double __expected[5][7] = {
                    {0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0},
                    {0, -1, -1, -1, -1, -1, -1},
                    {2, 2, 2, 2, 2, 2, 2},
                    {-1, -1, -1, -1, -1, -1, 0},
                };
                double* expected = flattened_row_major(5, 7, __expected);

                ASSERT(mat_equals(mat, expected, 5, 7));
                printf("Test Creation poisson matrix passed\n");
        }

        return 0;
}