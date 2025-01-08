/**
 * @file lib_poisson1D.c
 * @brief Numerical library developed to manipulate matrices and vectors for the Poisson 1D problem
 */

#include "lib_poisson1D.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv) {
        int filler = *kv;
        int nb_cols = *la;
        int nb_rows = *lab;

        for (int i = 0; i < nb_cols * nb_rows; i += nb_rows) {
                // Fill the free space
                for (int j = 0; j < filler; j++) {
                        AB[i + j] = 0;
                }

                // Each line has -1, 2, -1
                AB[i + filler] = -1;
                AB[i + filler + 1] = 2;
                AB[i + filler + 2] = -1;
        }

        // Make the final elem 0
        AB[nb_cols * nb_rows - 1] = 0;
        // Make the first elem of the first non kv line 0
        AB[filler] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv) {
        int filler = *kv;
        int nb_cols = *la;
        int nb_rows = *lab;

        // Fill everything with 0's
        for (int i = 0; i < nb_cols * nb_rows; i++) {
                AB[i] = 0;
        }

        // Fill the first last line with 1's (the diagonal)
        for (int i = 0; i < nb_cols; i++) {
                int idx = i * nb_rows + filler;
                AB[idx] = 1;
        }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
        int nb_points = *la;
        double T0 = *BC0;
        double T1 = *BC1;

        // Set the Dirichlet boundary conditions
        RHS[0] = T0;
        RHS[nb_points - 1] = T1;

        // Set the rest to 0
        for (int i = 1; i < nb_points - 1; i++) {
                RHS[i] = 0.0;
        }
}

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) {
        double T0 = *BC0;
        double T1 = *BC1;
        int nb_points = *la;

        for (int i = 0; i < nb_points; i++) {
                EX_SOL[i] = T0 + X[i] * (T1 - T0);
        }
}

void set_grid_points_1D(double* x, int* la) {
        int nb_points = *la;

        // Uniform grid
        double step = 1.0 / (1.0 * (nb_points + 1));
        for (int i = 0; i < nb_points; i++) {
                x[i] = (i + 1) * step;
        }
}

double relative_forward_error(double* x_exact, double* x_approx, int* la) {
        int nb_points = *la;
        double norm_x = sqrt(cblas_ddot(nb_points, x_exact, 1, x_exact, 1));   // ||x||
        cblas_daxpy(nb_points, -1, x_approx, 1, x_exact, 1);                   // x_exact = x_exact - x_approx
        double norm_res = sqrt(cblas_ddot(nb_points, x_exact, 1, x_exact, 1)); // ||x - x_approx||
        return norm_res / norm_x;                                              // ||x - x_approx|| / ||x||
}

// i and j ranges 0 to *la -1
int indexABCol(int i, int j, int* lab) {
        return (j + 1) * (*lab - 1) + i - 1;
}

int dgbtrftridiag(int* la, int* n, int* kl, int* ku, double* AB, int* lab, int* ipiv, int* info) {
        int nb_rows = *lab;
        int nb_cols = *la;
        int upper_diags = *ku;
        int lower_diags = *kl;
        int result_matrix_size = *n;

        // Make sure it's a tridiagonal matrix and it has room for the LU factorisation
        if (nb_rows != 4 || lower_diags != 1 || upper_diags != 1 || nb_cols < result_matrix_size) {
                *info = -1;
                return *info;
        }

        // Gaussian elimination with partial pivoting
        for (int k = 0; k < result_matrix_size - 1; k++) {
                // Find the biggest pivot
                // We optimise the search by only looking at the next row since it's a tridiagonal matrix
                // and looking further would only lead to 0's
                int pivot_idx;
                double pivot_value;
                if (fabs(AB[indexABCol(k, k, lab)]) > fabs(AB[indexABCol(k + 1, k, lab)])) {
                        pivot_idx = k;
                        pivot_value = AB[indexABCol(k, k, lab)];
                }
                else {
                        pivot_idx = k + 1;
                        pivot_value = AB[indexABCol(k + 1, k, lab)];
                }

                // Swap if necessary
                if (pivot_idx != k) {
                        // Swap rows k and pivot_idx
                        for (int j = 0; j < *lab; j++) {
                                double tmp = AB[indexABCol(pivot_idx, j, lab)];
                                AB[indexABCol(pivot_idx, j, lab)] = AB[indexABCol(k, j, lab)];
                                AB[indexABCol(k, j, lab)] = tmp;
                        }
                        // Mark the pivot
                        ipiv[k] = pivot_idx;
                }
                else {
                        // No row swap
                        ipiv[k] = k + 1;
                }

                // Perform the elimination, only the next row is affected since it's a tridiagonal matrix
                AB[indexABCol(k + 1, k, lab)] /= AB[indexABCol(k, k, lab)];
                AB[indexABCol(k + 1, k + 1, lab)] -= AB[indexABCol(k + 1, k, lab)] * AB[indexABCol(k, k + 1, lab)];
        }
        // Mark the last pivot
        ipiv[result_matrix_size - 1] = result_matrix_size;

        // Everything went well
        *info = 0;
        return *info;
}

void dcsrmv(char is_M_transposed, double alpha, CSRMatrix* M, double* x, size_t incx, double beta, double* y, size_t incy) {
        // If M is not transposed, we can directly do the dot product of each row of M with x and not have to use extra memory
        if (is_M_transposed == 'N') {
                for (int i = 0; i < M->nb_rows; i++) {
                        double sum = 0;
                        for (int nz_idx = M->row_ptr[i]; nz_idx < M->row_ptr[i + 1]; nz_idx++) {
                                sum += M->values[nz_idx] * x[M->col_index[nz_idx] * incx];
                        }
                        y[i * incy] = alpha * sum + beta * y[i * incy];
                }
        }
        // If M is transposed, we can't do it line by line, we need to store each result in their appropriate column and then sum them up
        else {
                double* result = (double*)calloc(M->nb_cols, sizeof(double));
                for (int i = 0; i < M->nb_rows; i++) {
                        for (int nz_idx = M->row_ptr[i]; nz_idx < M->row_ptr[i + 1]; nz_idx++) {
                                result[M->col_index[nz_idx]] += M->values[nz_idx] * x[i * incx];
                        }
                }
                for (int i = 0; i < M->nb_cols; i++) {
                        y[i * incy] = alpha * result[i] + beta * y[i * incy];
                }
                free(result);
        }
}

void dcscmv(char is_M_transposed, double alpha, CSCMatrix* M, double* x, size_t incx, double beta, double* y, size_t incy) {
        // If M is not transposed, we have to do it column by column, and then update the result
        if (is_M_transposed == 'N') {
                double* result = (double*)calloc(M->nb_rows, sizeof(double));
                for (int i = 0; i < M->nb_cols; i++) {
                        for (int nz_idx = M->col_ptr[i]; nz_idx < M->col_ptr[i + 1]; nz_idx++) {
                                result[M->row_index[nz_idx]] += M->values[nz_idx] * x[i * incx];
                        }
                }
                for (int i = 0; i < M->nb_rows; i++) {
                        y[i * incy] = alpha * result[i] + beta * y[i * incy];
                }
                free(result);
        }
        // If M is transposed, we can do it line by line, and directly update the result
        else {
                for (int i = 0; i < M->nb_cols; i++) {
                        double sum = 0;
                        for (int nz_idx = M->col_ptr[i]; nz_idx < M->col_ptr[i + 1]; nz_idx++) {
                                sum += M->values[nz_idx] * x[M->row_index[nz_idx] * incx];
                        }
                        y[i * incy] = alpha * sum + beta * y[i * incy];
                }
        }
}

// Computes the LU factorisation of a CSR matrix
int dcsrtrf(CSRMatrix* A, int* ipiv, int* info) {
        // Gaussian elimination with partial pivoting
        for (int k = 0; k < A->nb_rows - 1; k++) {
                // Find the biggest pivot
                int pivot_idx;
                double pivot_value;

                if (fabs(csr_elem_at(A, k, k)) > fabs(csr_elem_at(A, k + 1, k))) {
                        pivot_idx = k;
                        pivot_value = csr_elem_at(A, k, k);
                }
                else {
                        pivot_idx = k + 1;
                        pivot_value = csr_elem_at(A, k + 1, k);
                }

                // Swap if necessary
                if (pivot_idx != k) {
                        // Swap rows k and pivot_idx
                        for (int j = 0; j < A->nb_non_zero; j++) {
                                double tmp = csr_elem_at(A, pivot_idx, j);
                                csr_write_at(A, pivot_idx, j, csr_elem_at(A, k, j));
                                csr_write_at(A, k, j, tmp);
                        }
                        // Mark the pivot
                        ipiv[k] = pivot_idx;
                }
                else {
                        // No row swap
                        ipiv[k] = k + 1;
                }

                // Perform the elimination
                double factor = csr_elem_at(A, k + 1, k) / csr_elem_at(A, k, k);
                for (int j = 0; j < A->nb_non_zero; j++) {
                        csr_write_at(A, k + 1, j, csr_elem_at(A, k + 1, j) - factor * csr_elem_at(A, k, j));
                }
        }
        // Mark the last pivot
        ipiv[A->nb_rows - 1] = A->nb_rows;

        // Everything went well
        *info = 0;
        return *info;
}

void dcsrtrs(CSRMatrix* A, int* ipiv, double* b, double* x) {
        // Forward substitution
        for (int i = 0; i < A->nb_rows; i++) {
                x[i] = b[ipiv[i] - 1];
                for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                        x[i] -= A->values[j] * x[A->col_index[j]];
                }
        }

        // Backward substitution
        for (int i = A->nb_rows - 1; i >= 0; i--) {
                for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                        x[i] -= A->values[j] * x[A->col_index[j]];
                }
                x[i] /= A->values[A->row_ptr[i]];
        }
}

// Solve the linear system Ax = b by doing LU then forward and backward substitution
void dcsrsv(CSRMatrix* A, double* b, double* x) {
        int* ipiv = (int*)malloc(A->nb_rows * sizeof(int));
        int info;
        dcsrtrf(A, ipiv, &info);
        dcsrtrs(A, ipiv, b, x);

        // re-order x
        double* x_tmp = (double*)malloc(A->nb_rows * sizeof(double));
        for (int i = 0; i < A->nb_rows; i++) {
                x_tmp[ipiv[i] - 1] = x[i];
        }
        for (int i = 0; i < A->nb_rows; i++) {
                x[i] = x_tmp[i];
        }

        free(x_tmp);
        free(ipiv);
}