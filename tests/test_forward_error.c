/**
 * @file test_forward_error.c
 * @brief Test the relative forward error
 */

#include "atlas_headers.h"
#include "lib_poisson1D.h"
#include "test.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
        // Test de la forward error
        {
                int size = 5;
                double* x_exact = (double*)malloc(size * sizeof(double));
                double* x_approx = (double*)malloc(size * sizeof(double));

                for (int i = 0; i < size; i++) {
                        x_exact[i] = i;
                        x_approx[i] = i;
                }

                double forward_error = relative_forward_error(x_exact, x_approx, &size);
                ASSERT(forward_error < 1e-15);
                printf("Test forward error 1 passed\n");
                free(x_exact);
                free(x_approx);
        }

        {
                int size = 2;
                double* x_exact = (double*)malloc(size * sizeof(double));
                double* x_approx = (double*)malloc(size * sizeof(double));
                // x_exact = [1, 0], ||x_exact|| = 1
                x_exact[0] = 1;
                x_exact[1] = 0;

                // x_approx = [2, 1], x_approx - x = [1, 1], ||x_approx - x|| = sqrt(2)
                x_approx[0] = 2;
                x_approx[1] = 1;

                double forward_error = relative_forward_error(x_exact, x_approx, &size);
                ASSERT(fabs(forward_error - sqrt(2)) < 1e-6);
                printf("Test forward error 2 passed\n");
                free(x_exact);
                free(x_approx);
        }
        return 0;
}