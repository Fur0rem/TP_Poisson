/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include "atlas_headers.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);
double relative_forward_error(double* x, double* y, int* la);
void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename);
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void plot_convergence_history(double* vec, int size, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);
void eig_poisson1D(double* eigval, int* la);
double eigmax_poisson1D(int* la);
double eigmin_poisson1D(int* la);
double richardson_alpha_opt(int* la);
void richardson_alpha(double* AB, double* RHS, double* X, double* alpha_rich, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
		      double* resvec, int* nbite);
void extract_MB_jacobi_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv);
void extract_MB_gauss_seidel_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv);
void richardson_MB(double* AB, double* RHS, double* X, double* MB, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
		   double* resvec, int* nbite);
int indexABCol(int i, int j, int* lab);
int dgbtrftridiag(int* la, int* n, int* kl, int* ku, double* AB, int* lab, int* ipiv, int* info);

// CSR format
typedef struct CSRMatrix {
	int nb_rows;	 ///< Number of rows
	int nb_cols;	 ///< Number of columns
	int nb_non_zero; ///< Number of non zero elements
	double* values;	 ///< Values of the matrix
	int* col_index;	 ///< Column index of the non zero elements
	int* row_ptr;	 ///< Row pointer
} CSRMatrix;

double* csr_to_dense_col_major(CSRMatrix* csr);
CSRMatrix dense_col_major_to_csr(double* dense, int nb_lines, int nb_cols);

void print_csr_matrix(CSRMatrix* csr);

CSRMatrix poisson1D_csr_matrix(int nb_equations);

// Matrix in CSC format
typedef struct {
	int nb_rows;	 ///< Number of rows
	int nb_cols;	 ///< Number of columns
	int nb_non_zero; ///< Number of non zero elements
	double* values;	 ///< Values of the matrix
	int* row_index;	 ///< Row index of the non zero elements
	int* col_ptr;	 ///< Column pointer
} CSCMatrix;

double* csc_to_dense_col_major(CSCMatrix* csc);
CSCMatrix dense_col_major_to_csc(double* dense, int nb_lines, int nb_cols);

void print_csc_matrix(CSCMatrix* csc);

CSCMatrix poisson1D_csc_matrix(int nb_equations);

CSRMatrix csr_clone(CSRMatrix* csr);
CSCMatrix csc_clone(CSCMatrix* csc);

void csc_free(CSCMatrix* csc);
void csr_free(CSRMatrix* csr);

// alpha * M * x + beta * y
/**
 * @brief Compute the matrix-vector product y = alpha * M * x + beta * y
 *
 * @param[in] is_M_transposed 'N' if M is not transposed, 'T' if M is transposed
 * @param[in] alpha scalar alpha to multiply M * x by
 * @param[in] M matrix M in CSR format
 * @param[in] x vector x
 * @param[in] beta scalar beta to multiply y by
 * @param[in,out] y vector y, will be updated with the result
 * @param[in] incx increment for x
 * @param[in] incy increment for y
 */
void dcsrmv(char is_M_transposed, double alpha, CSRMatrix* M, double* x, double beta, double* y, size_t incx, size_t incy);