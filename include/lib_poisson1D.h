#ifndef LIB_POISSON1D_H
#define LIB_POISSON1D_H

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

/**
 * @brief Set the 1D Poisson matrix in GB format, column major
 * @param[out] AB Poisson matrix in GB format, assuming AB is allocated
 * @param[in] lab number of columns of AB
 * @param[in] la number of rows of AB
 * @param[in] kv number of superdiagonals of free space
 */
void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv);

/**
 * @brief Set the identity matrix in GB format, column major
 * @param[out] AB Identity matrix in GB format, assuming AB is allocated
 * @param[in] lab number of columns of AB
 * @param[in] la number of rows of AB
 * @param[in] kv number of superdiagonals of free space
 */
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv);

/**
 * @brief Set the system of equations right-hand side for the heat equation with Dirichlet boundary conditions
 * @param[out] RHS right-hand side of the system of equations, assuming RHS is allocated
 * @param[in] la number of equations in the problem to solve
 * @param[in] BC0 value of the Dirichlet boundary condition at the left boundary
 * @param[in] BC1 value of the Dirichlet boundary condition at the right boundary
 */
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);

/**
 * @brief Computes the analytical solution of the 1D Poisson problem with Dirichlet boundary conditions
 * @param[out] EX_SOL analytical solution of the 1D Poisson problem, assuming EX_SOL is allocated
 * @param[in] X grid points
 * @param[in] la number of equations in the problem to solve
 * @param[in] BC0 value of the Dirichlet boundary condition at the left boundary
 * @param[in] BC1 value of the Dirichlet boundary condition at the right boundary
 */
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);

/**
 * @brief Set the grid points for the 1D Poisson problem
 * @param[out] x grid points, assuming x is allocated
 * @param[in] la number of equations in the problem to solve
 */
void set_grid_points_1D(double* x, int* la);

/**
 * @brief Compute the relative forward error between two vectors
 * @param[in] x first vector
 * @param[in] y second vector
 * @param[in] la number of elements in the vectors
 * @return relative forward error
 */
double relative_forward_error(double* x, double* y, int* la);

/**
 * @brief Write a matrix in GB format to a file, with indices
 * @param[in] AB matrix in GB format
 * @param[in] la number of rows of AB
 * @param[in] lab number of columns of AB
 * @param[in] filename name of the file to write to
 */
void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename);

/**
 * @brief Write a matrix in GB format to a file, in row major
 * @param[in] AB matrix in GB format
 * @param[in] la number of rows of AB
 * @param[in] lab number of columns of AB
 * @param[in] filename name of the file to write to
 */
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename);

/**
 * @brief Write a matrix in GB format to a file, in column major
 * @param[in] AB matrix in GB format
 * @param[in] la number of rows of AB
 * @param[in] lab number of columns of AB
 * @param[in] filename name of the file to write to
 */
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);

/**
 * @brief Write a vector to a file
 * @param[in] vec vector to write
 * @param[in] la number of elements in the vector
 * @param[in] filename name of the file to write to
 */
void write_vec(double* vec, int* la, char* filename);

/**
 * @brief Plot the convergence history and saves it to a file
 * @param[in] vec convergence history to write
 * @param[in] size number of elements in the convergence history
 * @param[in] filename name of the file to write to
 */
void plot_convergence_history(double* vec, int size, char* filename);

/**
 * @brief Write two vectors to a file, a pair of elements per line
 * @param[in] vec first vector to write
 * @param[in] x second vector to write
 * @param[in] la number of elements in the vectors
 * @param[in] filename name of the file to write to
 */
void write_xy(double* vec, double* x, int* la, char* filename);

/**
 * @brief Computes the eigenvalues of the Poisson matrix
 * @param[put] eigval eigenvalues of the Poisson matrix
 * @param[in] la number of equations in the problem to solve
 */
void eig_poisson1D(double* eigval, int* la);

/**
 * @brief Computes the maximum eigenvalue of the Poisson matrix
 * @param[in] la number of equations in the problem to solve
 * @return maximum eigenvalue of the Poisson matrix
 */
double eigmax_poisson1D(int* la);

/**
 * @brief Computes the minimum eigenvalue of the Poisson matrix
 * @param[in] la number of equations in the problem to solve
 * @return minimum eigenvalue of the Poisson matrix
 */
double eigmin_poisson1D(int* la);

/**
 * @brief Computes the optimal value of alpha for the Richardson method
 * @param[in] la number of equations in the problem to solve
 * @return optimal value of alpha for the Richardson method
 */
double richardson_alpha_opt(int* la);

/**
 * @brief Solves Ax = b with the richardson method, with a scalar multiplier for iterating.
 * @param[in] AB matrix of the Poisson problem, in general band tri-diagonal form
 * @param[in] RHS right-hand side of the Poisson problem
 * @param[out] X solution of the Poisson problem
 * @param[in] alpha_rich value of alpha for the Richardson method
 * @param[in] lab number of rows of AB
 * @param[in] la number of columns of AB
 * @param[in] ku number of upper diagonals of AB
 * @param[in] kl number of lower diagonals of AB
 * @param[in] tol tolerance for the Richardson method
 * @param[in] maxit maximum number of iterations for the Richardson method
 * @param[out] resvec vector of the norm of residuals, already allocated of the size maxit
 * @param[out] nbite number of iterations taken
 */
void richardson_alpha_tridiag(double* AB, double* RHS, double* X, double* alpha_rich, int* lab, int* la, int* ku, int* kl, double* tol,
			      int* maxit, double* resvec, int* nbite);

/**
 * @brief Computes the iteration matrix for the Jacobi method in the Richardson algorithm, in general band tri-diagonal form
 * @param[in] AB matrix of the Poisson problem, in general band tri-diagonal form
 * @param[out] MB iteration matrix for the Jacobi method, in general band tri-diagonal form
 * @param[in] lab number of rows of AB
 * @param[in] la number of columns of AB
 * @param[in] ku number of upper diagonals of AB
 * @param[in] kl number of lower diagonals of AB
 * @param[in] kv number of diagonals of MB
 */
void extract_MB_jacobi_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv);

/**
 * @brief Computes the iteration matrix for the Gauss-Seidel method in the Richardson algorithm, in general band tri-diagonal form
 * @param[in] AB matrix of the Poisson problem, in general band tri-diagonal form
 * @param[out] MB iteration matrix for the Gauss-Seidel method, in general band tri-diagonal form
 * @param[in] lab number of rows of AB
 * @param[in] la number of columns of AB
 * @param[in] ku number of upper diagonals of AB
 * @param[in] kl number of lower diagonals of AB
 * @param[in] kv number of diagonals of MB
 */
void extract_MB_gauss_seidel_tridiag(double* AB, double* MB, int* lab, int* la, int* ku, int* kl, int* kv);

/**
 * @brief Solve the Poisson problem using the Richardson method with the iteration matrix MB
 * @param[in] AB matrix of the Poisson problem, in general band tri-diagonal form
 * @param[in] RHS right-hand side of the Poisson problem
 * @param[out] X solution of the Poisson problem, already allocated
 * @param[in] MB iteration matrix for the Richardson method, in general band tri-diagonal form
 * @param[in] lab number of rows of AB
 * @param[in] la number of columns of AB
 * @param[in] ku number of upper diagonals of AB
 * @param[in] kl number of lower diagonals of AB
 * @param[in] tol tolerance for the Richardson method
 * @param[in] maxit maximum number of iterations for the Richardson method
 * @param[out] resvec vector of the norm of residuals, already allocated of the size maxit
 * @param[out] nbite number of iterations
 */
void richardson_MB_tridiag(double* AB, double* RHS, double* X, double* MB, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit,
			   double* resvec, int* nbite);

/**
 * @brief Computes the index of the element in the AB matrix in column major format
 * @param[in] i row index
 * @param[in] j column index
 * @param[in] lab number of columns of AB
 * @return index of the element in the AB matrix in column major format
 */
int indexABCol(int i, int j, int* lab);

/**
 * @brief Compute the LU factorisation of a tridiagonal matrix

 */
int dgbtrftridiag(int* la, int* n, int* kl, int* ku, double* AB, int* lab, int* ipiv, int* info);

/**
 * @brief Compressed Sparse Row Matrix
 */
typedef struct CSRMatrix {
	int nb_rows;	 ///< Number of rows
	int nb_cols;	 ///< Number of columns
	int nb_non_zero; ///< Number of non zero elements
	double* values;	 ///< Values of the matrix
	int* col_index;	 ///< Column index of the non zero elements
	int* row_ptr;	 ///< Row pointer
} CSRMatrix;

/**
 * @brief Convert a CSR matrix to a dense matrix in column major format
 * @param[in] csr CSR matrix to convert
 * @return Dense matrix in column major format
 */
double* csr_to_dense_col_major(CSRMatrix* csr);

/**
 * @brief Convert a dense matrix in column major format to a CSR matrix
 * @param[in] dense Dense matrix in column major format
 * @param[in] nb_rows Number of lines of the dense matrix
 * @param[in] nb_cols Number of columns of the dense matrix
 * @return CSR matrix
 */
CSRMatrix dense_col_major_to_csr(double* dense, int nb_rows, int nb_cols);

/**
 * @brief Print a CSR matrix
 * @param[in] csr CSR matrix to print
 */
void print_csr_matrix(CSRMatrix* csr);

/**
 * @brief Create a 1D Poisson matrix in CSR format
 * @param[in] nb_equations Number of equations in the problem to solve
 * @return Poisson matrix in CSR format
 */
CSRMatrix poisson1D_csr_matrix(int nb_equations);

/**
 * @brief Compressed Sparse Column Matrix
 */
typedef struct {
	int nb_rows;	 ///< Number of rows
	int nb_cols;	 ///< Number of columns
	int nb_non_zero; ///< Number of non zero elements
	double* values;	 ///< Values of the matrix
	int* row_index;	 ///< Row index of the non zero elements
	int* col_ptr;	 ///< Column pointer
} CSCMatrix;

/**
 * @brief Convert a CSC matrix to a dense matrix in column major format
 * @param[in] csc CSC matrix to convert
 * @return Dense matrix in column major format
 */
double* csc_to_dense_col_major(CSCMatrix* csc);

/**
 * @brief Convert a dense matrix in column major format to a CSC matrix
 * @param[in] dense Dense matrix in column major format
 * @param[in] nb_rows Number of lines of the dense matrix
 * @param[in] nb_cols Number of columns of the dense matrix
 * @return CSC matrix
 */
CSCMatrix dense_col_major_to_csc(double* dense, int nb_rows, int nb_cols);

/**
 * @brief Print a CSC matrix
 * @param[in] csc CSC matrix to print
 */
void print_csc_matrix(CSCMatrix* csc);

/**
 * @brief Create a 1D Poisson matrix in CSC format
 * @param[in] nb_equations Number of equations in the problem to solve
 * @return Poisson matrix in CSC format
 */
CSCMatrix poisson1D_csc_matrix(int nb_equations);

/**
 * @brief Clone a CSR matrix
 * @param[in] csr CSR matrix to clone
 * @return Cloned CSR matrix
 */
CSRMatrix csr_clone(CSRMatrix* csr);

/**
 * @brief Clone a CSC matrix
 * @param[in] csc CSC matrix to clone
 * @return Cloned CSC matrix
 */
CSCMatrix csc_clone(CSCMatrix* csc);

/**
 * @brief Free a CSR matrix
 * @param[in] csr CSR matrix to free
 */
void csc_free(CSCMatrix* csc);

/**
 * @brief Free a CSC matrix
 * @param[in] csc CSC matrix to free
 */
void csr_free(CSRMatrix* csr);

/**
 * @brief Get the element at position (i, j) in a CSR matrix
 * @param[in] csr CSR matrix
 * @param[in] i Row index
 * @param[in] j Column index
 * @return Element at position (i, j)
 */
double csr_elem_at(CSRMatrix* csr, int i, int j);

/**
 * @brief Write a value at position (i, j) in a CSR matrix
 * @param[in] csr CSR matrix
 * @param[in] i Row index
 * @param[in] j Column index
 * @param[in] value Value to write
 */
void csr_write_at(CSRMatrix* csr, int i, int j, double value);

/**
 * @brief Get the element at position (i, j) in a CSC matrix
 * @param[in] csc CSC matrix
 * @param[in] i Row index
 * @param[in] j Column index
 * @return Element at position (i, j)
 */
double csc_elem_at(CSCMatrix* csc, int i, int j);

/**
 * @brief Write a value at position (i, j) in a CSC matrix
 * @param[in] csc CSC matrix
 * @param[in] i Row index
 * @param[in] j Column index
 * @param[in] value Value to write
 */
void csc_write_at(CSCMatrix* csc, int i, int j, double value);

/**
 * @brief Create a CSR matrix from a diagonal
 * @param[in] diag Diagonal values
 * @param[in] nb_elements Number of elements in the diagonal
 * @return CSR matrix
 */
CSRMatrix csr_from_diag(double* diag, int nb_elements);

/**
 * @brief Create a CSC matrix from a diagonal
 * @param[in] diag Diagonal values
 * @param[in] nb_elements Number of elements in the diagonal
 * @return CSC matrix
 */
CSCMatrix csc_from_diag(double* diag, int nb_elements);

/**
 * @brief Create a CSR matrix from a tridiagonal matrix
 * @param[in] diag Diagonal values
 * @param[in] lower Lower diagonal values
 * @param[in] upper Upper diagonal values
 * @param[in] nb_elements Number of elements in the diagonal
 * @return CSR matrix
 */
CSRMatrix csr_from_tridiag(double* diag, double* lower, double* upper, int nb_elements);

/**
 * @brief Create a CSC matrix from a tridiagonal matrix
 * @param[in] diag Diagonal values
 * @param[in] lower Lower diagonal values
 * @param[in] upper Upper diagonal values
 * @param[in] nb_elements Number of elements in the diagonal
 * @return CSC matrix
 */
CSCMatrix csc_from_tridiag(double* diag, double* lower, double* upper, int nb_elements);

/**
 * @brief Create a CSR matrix from a lower triangular matrix
 * @param[in] mat Lower triangular matrix
 * @param[in] nb_rows Number of rows
 * @return CSR matrix
 */
CSRMatrix csr_from_lower_triangular(double* mat, int nb_rows);

/**
 * @brief Create a CSC matrix from a lower triangular matrix
 * @param[in] mat Lower triangular matrix
 * @param[in] nb_rows Number of rows
 * @return CSC matrix
 */
CSCMatrix csc_from_lower_triangular(double* mat, int nb_rows);

/**
 * @brief Compute the matrix-vector product y = alpha * M * x + beta * y
 * @param[in] is_M_transposed 'N' if M is not transposed, 'T' if M is transposed
 * @param[in] alpha scalar alpha to multiply M * x by
 * @param[in] M matrix M in CSR format
 * @param[in] x vector x
 * @param[in] beta scalar beta to multiply y by
 * @param[in,out] y vector y, will be updated with the result
 * @param[in] incx increment for x
 * @param[in] incy increment for y
 */
void dcsrmv(char is_M_transposed, double alpha, CSRMatrix* M, double* x, size_t incx, double beta, double* y, size_t incy);

/**
 * @brief Extracts the Gauss-Seidel iteration matrix from the richardson algorithm in CSR format
 * @param[in] AB Equation matrix in CSR format
 * @return Gauss-Seidel iteration matrix in CSR format
 */
CSRMatrix extract_MB_gauss_seidel_csr(CSRMatrix* AB);

/**
 * @brief Extracts the Jacobi iteration matrix from the richardson algorithm in CSR format
 * @param[in] AB Equation matrix in CSR format
 * @return Jacobi iteration matrix in CSR format
 */
CSRMatrix extract_MB_jacobi_csr(CSRMatrix* AB);

/**
 * @brief Solves iteratively the system of linear equations Ax = B using the Richardson method with a step alpha
 * @param[in] AB Equation matrix in CSR format (A)
 * @param[in] RHS Right-hand side vector (B)
 * @param[out] X Solution vector (x)
 * @param[out] alpha_rich Optimal alpha for the Richardson method
 * @param[in] tol Tolerance for the residual
 * @param[in] maxit Maximum number of iterations
 * @param[out] resvec Residual vector for each iteration
 * @param[out] nbite Number of iterations taken
 */
void richardson_alpha_csr(CSRMatrix* AB, double* RHS, double* X, double* alpha_rich, double* tol, int* maxit, double* resvec, int* nbite);

/**
 * @brief Solves iteratively the system of linear equations Ax = B using the Richardson method with an iteration matrix
 * @param[in] AB Equation matrix in CSR format (A)
 * @param[in] RHS Right-hand side vector (B)
 * @param[out] X Solution vector (x)
 * @param[in] MB Iteration matrix in CSR format (M⁻¹)
 * @param[in] tol Tolerance for the residual
 * @param[in] maxit Maximum number of iterations
 * @param[out] resvec Residual vector for each iteration
 * @param[out] nbite Number of iterations taken
 */
void richardson_MB_csr(CSRMatrix* AB, double* RHS, double* X, CSRMatrix* MB, double* tol, int* maxit, double* resvec, int* nbite);

/**
 * @brief Triangularisation of a CSR matrix
 * @param[in] A Matrix in CSR format
 * @param[out] ipiv Pivot indices
 * @param[out] info Return status
 */
int dcsrtrf(CSRMatrix* A, int* ipiv, int* info);

/**
 * @brief Solve a system of linear equations Ax = b using the LU decomposition, assuming A has been triangularised using dcsrtrf
 * @param[in] A Matrix in CSR format, triangularised using dcsrtrf
 * @param[in] ipiv Pivot indices
 * @param[in,out] b Right-hand side vector
 * @param[out] x Solution vector
 */
void dcsrtrs(CSRMatrix* A, int* ipiv, double* b, double* x);

/**
 * @brief Solve a system of linear equations Ax = b using the LU decomposition
 * Calls dcsrtrf and dcsrtrs
 * @param[in] A Matrix in CSR format
 * @param[in,out] b Right-hand side vector
 * @param[out] x Solution vector
 */
void dcsrsv(CSRMatrix* A, double* b, double* x);

#endif // LIB_POISSON1D_H