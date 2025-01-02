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

// λ^k = 4 sin²(kπh / 2)
// eigval (out) : vector of eigenvalues to compute
// la (in) : size of the vector
void eig_poisson1D(double* eigval, int* la) {
	double h = 1.0 / (*la + 1);
	for (int k = 1; k <= *la; k++) {
		double sin_theta = sin(k * M_PI * h / 2.0);
		eigval[k - 1] = 4.0 * sin_theta * sin_theta;
	}
}

// compute_eigen_values : compute the eigenvalues of the Poisson matrix
double* compute_eigen_values(int* la) {
	double* eigval = (double*)malloc(*la * sizeof(double));
	eig_poisson1D(eigval, la);
	return eigval;
}

// eigmax_poisson1D : compute the maximum eigenvalue of the Poisson matrix
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

// eigmin_poisson1D : compute the minimum eigenvalue of the Poisson matrix
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

// richardson_alpha_opt : compute the optimal value of alpha for the Richardson method
// α_opt = 2 / (λmin + λmax)
double richardson_alpha_opt(int* la) {
	return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

double* residual(double* A, double* x, double* b, int nb_cols, int nb_lines) {
	double* res = (double*)malloc(nb_cols * sizeof(double));
	for (int i = 0; i < nb_cols; i++) {
		res[i] = b[i];
		for (int j = 0; j < nb_lines; j++) {
			res[i] -= A[i * nb_lines + j] * x[j];
		}
	}
	return res;
}

// richardson_alpha : solve the Poisson problem using the Richardson method
// AB (in, matrix) : the matrix of the Poisson problem, in general band tri diagonal form
// RHS (in, vec) : the right-hand side of the Poisson problem
// X (out, vec) : the solution of the Poisson problem, already allocated
// alpha_rich (in, scalar) : the value of alpha for the Richardson method
// lab (in, scalar) : number of lines of AB
// la (in, scalar) : number of columns of AB
// ku (in, scalar) : number of upper diagonals of AB
// kl (in, scalar) : number of lower diagonals of AB
// tol (in, scalar) : tolerance for the Richardson method
// maxit (in, scalar) : maximum number of iterations for the Richardson method
// resvec (out, vec) : vector of the norm of residuals, already allocated of the size maxit
// nbite (out, scalar) : number of iterations
void richardson_alpha(double* AB, double* RHS, double* X, double* alpha_rich, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
		      double* resvec, int* nbite) {
	double alpha_opt = *alpha_rich;
	int nb_lines = *lab;
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
	cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_lines, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter = 1;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + alpha * tmp_vec
		cblas_daxpy(nb_cols, alpha_opt, tmp_vec, 1, X, 1);

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_lines, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}
	}

	*nb_iters_final = iter;
	// residual_norms[iter] = norm_res;

	printf("Number of iterations: %d\n", iter);
	printf("X: ");
	for (int i = 0; i < nb_cols; i++) {
		printf("%lf ", X[i]);
	}

	free(tmp_vec);
}

// Extract the M iteration matrix for the Jacobi method
// AB (in, matrix) : the matrix of the Poisson problem, in general band tri diagonal form
// MB (out, matrix) : the iteration matrix for the Jacobi method, in general band tri diagonal form
// lab (in, scalar) : number of lines of AB
// la (in, scalar) : number of columns of AB
// ku (in, scalar) : number of upper diagonals of AB
// kl (in, scalar) : number of lower diagonals of AB
// kv (in, scalar) : number of diagonals of MB

// M-1 = D^-1 (diagonal)
void extract_MB_jacobi_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv) {
	int nb_lines = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int diagonals = *kv;

	// Set MB to 0
	memset(MB, 0, nb_lines * diagonals * sizeof(double));

	// Copy the diagonal of AB to MB
	for (int i = 0; i < nb_cols; i++) {
		MB[1 + i * (diagonals + upper_diags + lower_diags)] = 1.0 / AB[1 + i * (diagonals + upper_diags + lower_diags)];
	}
}

// Extract the M iteration matrix for the Gauss-Seidel method
// AB (in, matrix) : the matrix of the Poisson problem, in general band tri diagonal form
// MB (out, matrix) : the iteration matrix for the Gauss-Seidel method, in general band tri diagonal form
// lab (in, scalar) : number of lines of AB
// la (in, scalar) : number of columns of AB
// ku (in, scalar) : number of upper diagonals of AB
// kl (in, scalar) : number of lower diagonals of AB
// kv (in, scalar) : number of diagonals of MB

// M = D (diagonal) - E (lower triangle)
void extract_MB_gauss_seidel_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv) {
	int nb_lines = *lab;
	int nb_cols = *la;
	int upper_diags = *ku;
	int lower_diags = *kl;
	int diagonals = *kv;

	printf("diagonals : %d\n", diagonals);

	// Set MB to 0
	memset(MB, 0, nb_lines * diagonals * sizeof(double));

	// FIXME: arbitrary values for now
	// Diagonal is equal to 0.5
	for (int i = 0; i < nb_cols; i++) {
		MB[1 + i * (diagonals + upper_diags + lower_diags)] = 1.0 / AB[1 + i * (diagonals + upper_diags + lower_diags)];
	}

	// Lower diagonal is equal to 0.25
	for (int i = 1; i < nb_cols; i++) {
		MB[0 + i * (diagonals + upper_diags + lower_diags)] = 0.25;
	}

	printf("MB : ");
	for (int i = 0; i < nb_cols; i++) {
		printf("%lf ", MB[i]);
	}
	printf("\n");
}

// richardson_MB : solve the Poisson problem using the Richardson method with the iteration matrix MB
// AB (in, matrix) : the matrix of the Poisson problem, in general band tri diagonal form
// RHS (in, vec) : the right-hand side of the Poisson problem
// X (out, vec) : the solution of the Poisson problem, already allocated
// MB (in, matrix) : the iteration matrix for the Richardson method, in general band tri diagonal form, assumed it is already inverted
// lab (in, scalar) : number of lines of AB
// la (in, scalar) : number of columns of AB
// ku (in, scalar) : number of upper diagonals of AB
// kl (in, scalar) : number of lower diagonals of AB
// tol (in, scalar) : tolerance for the Richardson method
// maxit (in, scalar) : maximum number of iterations for the Richardson method
// resvec (out, vec) : vector of the norm of residuals, already allocated of the size maxit
// nbite (out, scalar) : number of iterations

// Used to compute the jacobi method or gauss-seidel method depending on the MB matrix
void richardson_MB(double* AB, double* RHS, double* X, double* MB, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
		   double* resvec, int* nbite) {
	// double alpha_opt = *alpha_rich;
	int nb_lines = *lab;
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
	cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_lines, X, 1, 1.0, tmp_vec, 1);

	// Compute the initial residual norm
	double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
	residual_norms[0] = norm_res;

	// Iterate until you reach the maximum number of iterations and break early if tolerance is reached
	int iter = 1;
	for (iter = 1; iter < max_iters; ++iter) {
		// X = X + alpha * tmp_vec
		// cblas_daxpy(nb_cols, alpha_opt, tmp_vec, 1, X, 1);

		// X = X + MB * tmp_vec
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, 1.0, MB, nb_lines, tmp_vec, 1, 1.0, X, 1);

		// printf("\nX k : ");
		// for (int i = 0; i < nb_cols; i++) {
		// 	printf("%lf ", X[i]);
		// }
		// printf("\n");

		// tmp_vec = RHS
		memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

		// tmp_vec = tmp_vec - AB * X
		cblas_dgbmv(
		    CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_lines, X, 1, 1.0, tmp_vec, 1);

		// Compute the residual and store it
		norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
		residual_norms[iter] = norm_res;

		// Check for convergence
		if (norm_res < tolerance) {
			break;
		}

		// printf("residual at iter %d is %lf\n", iter, residual_norms[iter]);
		// printf("X k+1 : ");
		// for (int i = 0; i < nb_cols; i++) {
		// 	printf("%lf ", X[i]);
		// }
		// printf("\n\n");
	}

	*nb_iters_final = iter;
	// residual_norms[iter] = norm_res;

	printf("Number of iterations: %d\n", iter);
	printf("X: ");
	for (int i = 0; i < nb_cols; i++) {
		printf("%lf ", X[i]);
	}

	free(tmp_vec);
}