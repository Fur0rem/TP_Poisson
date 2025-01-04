#include "atlas_headers.h"
#include "lib_poisson1D.h"
#include "test.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
	printf("--------- Test format CSR et CSC ---------\n");

	{
		CSRMatrix poisson = poisson1D_csr_matrix(4);
		print_csr_matrix(&poisson);
		double* dense = csr_to_dense_col_major(&poisson);
		print_mat_col_major(dense, poisson.nb_rows, poisson.nb_cols);
		double expected[4][4] = {
		    {2, -1, 0, 0},
		    {-1, 2, -1, 0},
		    {0, -1, 2, -1},
		    {0, 0, -1, 2},
		};
		double* expected_dense = flattened_row_major(4, 4, expected);
		ASSERT(mat_equals(dense, expected_dense, 4, 4));
		printf("Test format CSR passed\n");
	}
	{
		CSCMatrix poisson = poisson1D_csc_matrix(4);
		print_csc_matrix(&poisson);
		double* dense = csc_to_dense_col_major(&poisson);
		print_mat_col_major(dense, poisson.nb_rows, poisson.nb_cols);
		double expected[4][4] = {
		    {2, -1, 0, 0},
		    {-1, 2, -1, 0},
		    {0, -1, 2, -1},
		    {0, 0, -1, 2},
		};
		double* expected_dense = flattened_row_major(4, 4, expected);
		ASSERT(mat_equals(dense, expected_dense, 4, 4));
		printf("Test format CSC passed\n");
	}

	// CSR Tests
	{
		// Example taken from https://mathinsight.org/matrix_vector_multiplication
		double dense_mat[2][3] = {
		    {1, -1, 2},
		    {0, -3, 1},
		};
		double vec_x[6] = {2, 1, 0};
		double vec_y[2] = {0, 0};
		double* dense = flattened_row_major(2, 3, dense_mat);
		CSRMatrix csr = dense_col_major_to_csr(dense, 2, 3);
		dcsrmv('N', 1.0, &csr, vec_x, 0.0, vec_y, 1, 1);

		double expected[2] = {1, -3};
		ASSERT(vec_equals(vec_y, expected, 2));
	}

	{
		// Same example as above but with different alpha, beta, and increments
		double dense_mat[2][3] = {
		    {1, -1, 2},
		    {0, -3, 1},
		};
		double vec_x[6] = {2, 0, 1, 0, 0, 0};
		double vec_y[2] = {1, 1};
		CSRMatrix csr = dense_col_major_to_csr(flattened_row_major(2, 3, dense_mat), 2, 3);
		dcsrmv('N', 2.0, &csr, vec_x, 3.0, vec_y, 2, 1);
		double expected[2] = {5, -3};
		ASSERT(vec_equals(vec_y, expected, 2));
	}

	{
		// Same example as above but transposed
		double dense_mat[3][2] = {
		    {1, 0},
		    {-1, -3},
		    {2, 1},
		};
		double vec_x[6] = {2, 0, 1, 0, 0, 0};
		double vec_y[2] = {1, 1};
		CSRMatrix csr = dense_col_major_to_csr(flattened_row_major(3, 2, dense_mat), 3, 2);
		dcsrmv('T', 2.0, &csr, vec_x, 3.0, vec_y, 2, 1);
		double expected[2] = {5, -3};
		printf("vec_y = %lf %lf\n", vec_y[0], vec_y[1]);

		ASSERT(vec_equals(vec_y, expected, 2));
	}

	// CSC tests
	{
		// Example taken from https://mathinsight.org/matrix_vector_multiplication
		double dense_mat[2][3] = {
		    {1, -1, 2},
		    {0, -3, 1},
		};
		double vec_x[6] = {2, 1, 0};
		double vec_y[2] = {0, 0};
		double* dense = flattened_row_major(2, 3, dense_mat);
		CSCMatrix csc = dense_col_major_to_csc(dense, 2, 3);
		dcscmv('N', 1.0, &csc, vec_x, 0.0, vec_y, 1, 1);

		double expected[2] = {1, -3};
		ASSERT(vec_equals(vec_y, expected, 2));
	}

	{
		// Same example as above but with different alpha, beta, and increments
		double dense_mat[2][3] = {
		    {1, -1, 2},
		    {0, -3, 1},
		};
		double vec_x[6] = {2, 0, 1, 0, 0, 0};
		double vec_y[2] = {1, 1};
		CSCMatrix csc = dense_col_major_to_csc(flattened_row_major(2, 3, dense_mat), 2, 3);
		dcscmv('N', 2.0, &csc, vec_x, 3.0, vec_y, 2, 1);
		double expected[2] = {5, -3};
		ASSERT(vec_equals(vec_y, expected, 2));
	}

	{
		// Same example as above but transposed
		double dense_mat[3][2] = {
		    {1, 0},
		    {-1, -3},
		    {2, 1},
		};
		double vec_x[6] = {2, 0, 1, 0, 0, 0};
		double vec_y[2] = {1, 1};
		CSCMatrix csc = dense_col_major_to_csc(flattened_row_major(3, 2, dense_mat), 3, 2);
		dcscmv('T', 2.0, &csc, vec_x, 3.0, vec_y, 2, 1);
		double expected[2] = {5, -3};
		printf("vec_y = %lf %lf\n", vec_y[0], vec_y[1]);

		ASSERT(vec_equals(vec_y, expected, 2));
	}
}