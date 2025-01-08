/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>

#define ALPHA 0
#define JAC   1
#define GS    2

#define TRIDIAG 0
#define CSR	1
#define CSC	2

int main(int argc, char* argv[]) {
	printf("--------- ITERATIVE METHODS ---------\n\n");
	int ierr;
	int jj;
	int nbpoints, la;
	int ku, kl, lab, kv;
	int* ipiv;
	int info;
	int NRHS;
	int IMPLEM = 0;
	int FORMAT = TRIDIAG;
	double T0, T1;
	double *RHS, *SOL, *EX_SOL, *X;
	double* AB;
	double* MB;

	double temp, relres;

	double opt_alpha;

	if (argc == 3) {
		IMPLEM = atoi(argv[1]);
		FORMAT = atoi(argv[2]);
		printf("IMPLEM = %d, FORMAT = %d\n", IMPLEM, FORMAT);
	}
	else {
		// TODO: better error message
		perror("Usage: ./tp2_poisson1D_iter <IMPLEM> <FORMAT>, format=0 -> tridiag, format=1 -> csr, format=2 -> csc");
		exit(1);
	}

	/* Size of the problem */
	NRHS = 1;
	nbpoints = 12;
	la = nbpoints - 2;

	/* Dirichlet Boundary conditions */
	T0 = 5.0;
	T1 = 20.0;

	printf("--------- Poisson 1D ---------\n\n");
	RHS = (double*)malloc(sizeof(double) * la);
	SOL = (double*)calloc(la, sizeof(double));
	EX_SOL = (double*)malloc(sizeof(double) * la);
	X = (double*)malloc(sizeof(double) * la);

	/* Setup the Poisson 1D problem */
	/* General Band Storage */
	set_grid_points_1D(X, &la);
	set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
	set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

	write_vec(RHS, &la, "RHS.dat");
	write_vec(EX_SOL, &la, "EX_SOL.dat");
	write_vec(X, &la, "X_grid.dat");

	kv = 0;
	ku = 1;
	kl = 1;
	lab = kv + kl + ku + 1;

	AB = (double*)malloc(sizeof(double) * lab * la);
	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

	/* uncomment the following to check matrix A */
	write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

	/********************************************/
	/* Solution (Richardson with optimal alpha) */

	/* Computation of optimum alpha */

	/* Solve */
	double tol = 1e-3;
	int maxit = 1000;
	double* resvec;
	int nbite = 0;

	resvec = (double*)calloc(maxit, sizeof(double));

	/* Solve with Richardson alpha */
	if (FORMAT == TRIDIAG) {
		printf("Format : Tri-Diag\n");
		if (IMPLEM == ALPHA) {
			opt_alpha = richardson_alpha_opt(&la);
			printf("Optimal alpha for simple Richardson iteration is : %lf\n", opt_alpha);
			richardson_alpha_tridiag(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
			plot_convergence_history(resvec, nbite, "rapport/convergence_richardson_alpha_format_tri-diag");
			printf("Number of iterations for Richardson with optimal alpha is : %d\n", nbite);
		}
		/* Richardson General Tridiag */

		/* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
		kv = 1;
		ku = 1;
		kl = 1;
		MB = (double*)malloc(sizeof(double) * (lab)*la);
		if (IMPLEM == JAC) {
			printf("Jacobi\n");
			extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
		}
		else if (IMPLEM == GS) {
			printf("Gauss-Seidel\n");
			extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
		}

		/* Solve with General Richardson */
		if (IMPLEM == JAC || IMPLEM == GS) {
			write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
			richardson_MB_tridiag(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
			if (IMPLEM == JAC) {
				plot_convergence_history(resvec, nbite, "rapport/convergence_jacobi_format_tri-diag");
				printf("Number of iterations for Jacobi is : %d\n", nbite);
			}
			else if (IMPLEM == GS) {
				plot_convergence_history(resvec, nbite, "rapport/convergence_gauss-seidel_format_tri-diag");
				printf("Number of iterations for Gauss-Seidel is : %d\n", nbite);
			}
		}
	}
	else if (FORMAT == CSR) {
		printf("Format : CSR\n");
		// Extract CSR matrix
		CSRMatrix AB = poisson1D_csr_matrix(la);

		// printf("AB matrix:\n");
		// print_csr_matrix(&AB);

		// Solve with Richardson alpha
		if (IMPLEM == ALPHA) {
			opt_alpha = richardson_alpha_opt(&la);
			printf("Optimal alpha for simple Richardson iteration is : %lf\n", opt_alpha);
			richardson_alpha_csr(&AB, RHS, SOL, &opt_alpha, &tol, &maxit, resvec, &nbite);
			plot_convergence_history(resvec, nbite, "rapport/convergence_richardson_alpha_format_CSR");
			printf("Number of iterations for Richardson with optimal alpha is : %d\n", nbite);
		}
		/* Richardson General CSR */

		/* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
		CSRMatrix MB;
		if (IMPLEM == JAC) {
			printf("Jacobi CSR\n");
			MB = extract_MB_jacobi_csr(&AB);
		}
		else if (IMPLEM == GS) {
			printf("Gauss-Seidel CSR\n");
			MB = extract_MB_gauss_seidel_csr(&AB);
		}

		// print_csr_matrix(&MB);

		/* Solve with General Richardson */
		if (IMPLEM == JAC || IMPLEM == GS) {
			richardson_MB_csr(&AB, RHS, SOL, &MB, &tol, &maxit, resvec, &nbite);
			if (IMPLEM == JAC) {
				plot_convergence_history(resvec, nbite, "rapport/convergence_jacobi_format_CSR");
				printf("Number of iterations for Jacobi is : %d\n", nbite);
			}
			else if (IMPLEM == GS) {
				plot_convergence_history(resvec, nbite, "rapport/convergence_gauss-seidel_format_CSR");
				printf("Number of iterations for Gauss-Seidel is : %d\n", nbite);
			}
		}
	}
	else if (IMPLEM == CSC) {
		printf("Not implemented yet\n");
		exit(EXIT_FAILURE);
	}

	// Print x and exact solution
	printf("X found:\t");
	for (int i = 0; i < la; i++) {
		printf("%lf\t", SOL[i]);
	}
	printf("\n");

	printf("Exact exact:\t");
	for (int i = 0; i < la; i++) {
		printf("%lf\t", EX_SOL[i]);
	}
	printf("\n");

	/* Write solution */
	write_vec(SOL, &la, "SOL.dat");

	/* Write convergence history */
	write_vec(resvec, &nbite, "RESVEC.dat");

	free(resvec);
	free(RHS);
	free(SOL);
	free(EX_SOL);
	free(X);
	free(AB);
	free(MB);
	printf("\n\n--------- End -----------\n");
}
