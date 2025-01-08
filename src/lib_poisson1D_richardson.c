/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void eig_poisson1D(double* eigval, int* la) {
	double h = 1.0 / (*la + 1);
	// λ^k = 4 sin²(kπh / 2)
	for (int k = 1; k <= *la; k++) {
		double sin_theta = sin(k * M_PI * h / 2.0);
		eigval[k - 1] = 4.0 * sin_theta * sin_theta;
	}
}

double* compute_eigen_values(int* la) {
	double* eigval = (double*)malloc(*la * sizeof(double));
	eig_poisson1D(eigval, la);
	return eigval;
}

double eigmax_poisson1D(int* la) {
	double* eigval = compute_eigen_values(la);

	// Find the maximum eigenvalue
	double eigmax = eigval[0];
	for (int i = 1; i < *la; i++) {
		if (eigval[i] > eigmax) {
			eigmax = eigval[i];
		}
	}

	free(eigval);
	return eigmax;
}

double eigmin_poisson1D(int* la) {
	double* eigval = compute_eigen_values(la);

	// Find the minimum eigenvalue
	double eigmin = eigval[0];
	for (int i = 1; i < *la; i++) {
		if (eigval[i] < eigmin) {
			eigmin = eigval[i];
		}
	}

	free(eigval);
	return eigmin;
}

double richardson_alpha_opt(int* la) {
	// α_opt = 2 / (λmin + λmax)
	return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

double* residual(double* A, double* x, double* b, int nb_cols, int nb_rows) {
	double* res = (double*)malloc(nb_cols * sizeof(double));
	for (int i = 0; i < nb_cols; i++) {
		res[i] = b[i];
		// Difference between actual result and computed result
		for (int j = 0; j < nb_rows; j++) {
			res[i] -= A[i * nb_rows + j] * x[j];
		}
	}
	return res;
}

void richardson_alpha_tridiag(double* AB, double* RHS, double* X, double* alpha_rich, int* lab, int* la, int* ku, int* kl, double* tol,
			      int* maxit, double* resvec, int* nbite) {
	double alpha = *alpha_rich;
	int nb_rows = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int max_iters = *maxit;
	double tolerance = *tol;
	int* nb_iters_final = nbite;
	double* residual_norms = resvec;

	double* tmp_vec = malloc(sizeof(double) * nb_cols);
	if (tmp_vec == NULL) {
		fprintf(stderr, "Memory allocation failed for temp\n");
		exit(EXIT_FAILURE);
	}

	// tmp_vec = RHS
	memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

	// tmp_vec = tmp_vec - AB * X
	cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + alpha * tmp_vec
		cblas_daxpy(nb_cols, alpha, tmp_vec, 1, X, 1);

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}
	}

	*nb_iters_final = iter;

	free(tmp_vec);
}

void richardson_alpha_csr(CSRMatrix* AB, double* RHS, double* X, double* alpha_rich, double* tol, int* maxit, double* resvec, int* nbite) {
	double alpha = *alpha_rich;
	int nb_rows = AB->nb_rows;
	int nb_cols = AB->nb_cols;
	int max_iters = *maxit;
	double tolerance = *tol;
	int* nb_iters_final = nbite;
	double* residual_norms = resvec;

	double* tmp_vec = malloc(sizeof(double) * nb_cols);
	if (tmp_vec == NULL) {
		fprintf(stderr, "Memory allocation failed for temp\n");
		exit(EXIT_FAILURE);
	}

	// tmp_vec = RHS
	memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

	// tmp_vec = tmp_vec - AB * X
	dcsrmv('N', -1.0, AB, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + alpha * tmp_vec
		cblas_daxpy(nb_cols, alpha, tmp_vec, 1, X, 1);

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		dcsrmv('N', -1.0, AB, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}
	}

	*nb_iters_final = iter;

	free(tmp_vec);
}

// M^-1 = D^-1
void extract_MB_jacobi_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv) {
	int nb_rows = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int diagonals = *kv;

	// All elements are 0 except the diagonal
	memset(MB, 0, nb_rows * diagonals * sizeof(double));

	// Since D is diagonal, we can just invert the diagonal elements
	for (int i = 0; i < nb_cols; i++) {
		MB[1 + i * (diagonals + upper_diags + lower_diags)] = 1.0 / AB[1 + i * (diagonals + upper_diags + lower_diags)];
	}
}

// M^-1 = D^-1
CSRMatrix extract_MB_jacobi_csr(CSRMatrix* AB) {
	// Place for the diagonal elements
	CSRMatrix MB = {
	    .nb_rows = AB->nb_rows,
	    .nb_cols = AB->nb_cols,
	    .nb_non_zero = AB->nb_cols,
	    .values = (double*)malloc(AB->nb_cols * sizeof(double)),
	    .row_ptr = (int*)malloc((AB->nb_rows + 1) * sizeof(int)),
	    .col_index = (int*)malloc(AB->nb_cols * sizeof(int)),
	};
	MB.row_ptr[0] = 0;

	// Extract the diagonal of AB and invert it
	for (int i = 0; i < AB->nb_cols; i++) {
		MB.values[i] = 1.0 / csr_elem_at(AB, i, i);
		MB.col_index[i] = i;
		MB.row_ptr[i + 1] = i + 1;
	}

	return MB;
}

/**
 * @brief Compute the index of the element in the AB matrix
 * @param[in] i row index
 * @param[in] j column index
 * @param[in] lab number of columns of AB
 * @return index of the element in the AB matrix
 */
int idx_into_lower_triangular(int i, int j, int n) {
	return (i + 1) * i / 2 + j;
}

/**
 * @brief Invert a lower triangular matrix
 * @param[in] lower_triangular Lower triangular matrix
 * @param[out] inv_lower_triangular Inverted lower triangular matrix
 * @param[in] n Number of rows
 */
void invert_lower_triangular(double* lower_triangular, double* inv_lower_triangular, int n) {
	// Set the diagonal to reciprocal
	for (int i = 0; i < n; i++) {
		inv_lower_triangular[idx_into_lower_triangular(i, i, n)] = 1.0 / lower_triangular[idx_into_lower_triangular(i, i, n)];
	}

	// Now we iteratively compute the other elements from top to bottom such that A * inv_A = I
	for (int i = 1; i < n; i++) {
		// Compute the elements of the i-th row
		for (int j = 0; j < i; j++) {
			double sum = 0.0;
			for (int k = j; k < i; k++) {
				sum += lower_triangular[idx_into_lower_triangular(i, k, n)] *
				       inv_lower_triangular[idx_into_lower_triangular(k, j, n)];
			}
			inv_lower_triangular[idx_into_lower_triangular(i, j, n)] =
			    -sum / lower_triangular[idx_into_lower_triangular(i, i, n)];
		}
	}
}

// M = D (diagonal) - E (lower triangle)
// M^-1 = (D - E)^-1
void extract_MB_gauss_seidel_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv) {
	int nb_rows = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int diagonals = *kv;

	// Go from tri-diagonal to lower triangular (D - E)
	double* lower_triangular = (double*)malloc(((nb_cols + 1) * nb_cols / 2) * sizeof(double));
	memset(lower_triangular, 0, ((nb_cols + 1) * nb_cols / 2) * sizeof(double));

	// Copy the diagonal
	int current_idx = 0;
	for (int i = 0; i < nb_cols; i++) {
		lower_triangular[current_idx] = AB[1 + i * (diagonals + upper_diags + lower_diags)];
		current_idx += i + 2;
	}

	// Copy the lower diagonal
	for (int i = 1; i < nb_cols; i++) {
		lower_triangular[idx_into_lower_triangular(i, i - 1, nb_cols)] = AB[0 + i * (diagonals + upper_diags + lower_diags)];
	}

	// Invert the lower triangular matrix
	double* inv_lower_triangular = (double*)malloc(((nb_cols + 1) * nb_cols / 2) * sizeof(double));
	invert_lower_triangular(lower_triangular, inv_lower_triangular, nb_cols);

	// Truncate the inv_lower_triangular back into a tridiagonal matrix
	memset(MB, 0, nb_rows * (diagonals + upper_diags + lower_diags) * sizeof(double));
	current_idx = 0;

	// Copy the diagonal
	for (int i = 0; i < nb_cols; i++) {
		MB[1 + i * (diagonals + upper_diags + lower_diags)] = inv_lower_triangular[current_idx];
		current_idx += i + 2;
	}

	// Copy the lower diagonal
	for (int i = 1; i < nb_cols; i++) {
		MB[0 + i * (diagonals + upper_diags + lower_diags)] = inv_lower_triangular[idx_into_lower_triangular(i, i - 1, nb_cols)];
	}

	free(lower_triangular);
}

// M = D (diagonal) - E (lower triangle)
// M^-1 = (D - E)^-1
CSRMatrix extract_MB_gauss_seidel_csr(CSRMatrix* AB) {
	// Convert to dense then to lower triangular
	double* dense = csr_to_dense_col_major(AB);
	double* lower_triangular = (double*)malloc((AB->nb_rows + 1) * AB->nb_rows / 2 * sizeof(double));

	// Copy the diagonal
	int current_idx = 0;
	for (int i = 0; i < AB->nb_rows; i++) {
		lower_triangular[current_idx] = dense[i + i * AB->nb_rows];
		current_idx += i + 2;
	}

	// Copy the lower diagonal
	for (int i = 1; i < AB->nb_rows; i++) {
		lower_triangular[idx_into_lower_triangular(i, i - 1, AB->nb_rows)] = dense[i + (i - 1) * AB->nb_rows];
	}

	// Invert the lower triangular matrix
	double* inv_lower_triangular = (double*)malloc((AB->nb_rows + 1) * AB->nb_rows / 2 * sizeof(double));
	invert_lower_triangular(lower_triangular, inv_lower_triangular, AB->nb_rows);

	// Put it back into the CSR format
	CSRMatrix inv_lower_triangular_csr = csr_from_lower_triangular(inv_lower_triangular, AB->nb_rows);

	// Free memory
	free(dense);
	free(lower_triangular);
	free(inv_lower_triangular);

	return inv_lower_triangular_csr;
}

void richardson_MB_tridiag(double* AB, double* RHS, double* X, double* MB, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
			   double* resvec, int* nbite) {
	int nb_rows = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int max_iters = *maxit;
	double tolerance = *tol;
	int* nb_iters_final = nbite;
	double* residual_norms = resvec;

	// X = 0
	memset(X, 0, nb_cols * sizeof(double));

	double* tmp_vec = malloc(sizeof(double) * nb_cols);
	if (tmp_vec == NULL) {
		fprintf(stderr, "Memory allocation failed for temp\n");
		exit(EXIT_FAILURE);
	}

	// tmp_vec = RHS
	memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

	// tmp_vec = tmp_vec - AB * X
	cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + MB * tmp_vec
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, 1.0, MB, nb_rows, tmp_vec, 1, 1.0, X, 1);

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}
	}

	*nb_iters_final = iter;

	free(tmp_vec);
}

void richardson_MB_csr(CSRMatrix* AB, double* RHS, double* X, CSRMatrix* MB, double* tol, int* maxit, double* resvec, int* nbite) {
	int nb_rows = AB->nb_rows;
	int nb_cols = AB->nb_cols;
	int max_iters = *maxit;
	double tolerance = *tol;
	int* nb_iters_final = nbite;
	double* residual_norms = resvec;

	// X = 0
	memset(X, 0, nb_cols * sizeof(double));

	double* tmp_vec = malloc(sizeof(double) * nb_cols);
	if (tmp_vec == NULL) {
		fprintf(stderr, "Memory allocation failed for temp\n");
		exit(EXIT_FAILURE);
	}

	// tmp_vec = RHS
	memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

	// tmp_vec = tmp_vec - AB * X
	dcsrmv('N', -1.0, AB, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + MB * tmp_vec
		dcsrmv('N', 1.0, MB, tmp_vec, 1, 1.0, X, 1);

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		dcsrmv('N', -1.0, AB, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}
	}

	*nb_iters_final = iter;

	free(tmp_vec);
}