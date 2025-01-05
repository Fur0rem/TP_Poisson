// Matrix in CSR format
#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double* csr_to_dense_col_major(CSRMatrix* csr) {
	double* dense = (double*)malloc(csr->nb_rows * csr->nb_cols * sizeof(double));
	memset(dense, 0, csr->nb_rows * csr->nb_cols * sizeof(double));
	for (int i = 0; i < csr->nb_rows; i++) {
		for (int j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++) {
			dense[i + csr->col_index[j] * csr->nb_rows] = csr->values[j];
		}
	}
	return dense;
}

CSRMatrix dense_col_major_to_csr(double* dense, int nb_rows, int nb_cols) {
	CSRMatrix csr;
	csr.nb_rows = nb_rows;
	csr.nb_cols = nb_cols;
	csr.nb_non_zero = 0;
	for (int i = 0; i < nb_rows; i++) {
		for (int j = 0; j < nb_cols; j++) {
			if (dense[i + j * nb_rows] != 0) {
				csr.nb_non_zero++;
			}
		}
	}
	csr.values = (double*)malloc(csr.nb_non_zero * sizeof(double));
	csr.col_index = (int*)malloc(csr.nb_non_zero * sizeof(int));
	csr.row_ptr = (int*)malloc((nb_rows + 1) * sizeof(int));
	int non_zero_idx = 0;
	for (int i = 0; i < nb_rows; i++) {
		csr.row_ptr[i] = non_zero_idx;
		for (int j = 0; j < nb_cols; j++) {
			if (dense[i + j * nb_rows] != 0) {
				csr.values[non_zero_idx] = dense[i + j * nb_rows];
				csr.col_index[non_zero_idx] = j;
				non_zero_idx++;
			}
		}
	}
	csr.row_ptr[nb_rows] = non_zero_idx;
	return csr;
}

void print_csr_matrix(CSRMatrix* csr) {
	for (int i = 0; i < csr->nb_rows; i++) {
		for (int j = 0; j < csr->nb_cols; j++) {
			int found = 0;
			for (int k = csr->row_ptr[i]; k < csr->row_ptr[i + 1]; k++) {
				if (csr->col_index[k] == j) {
					printf("%.4lf\t", csr->values[k]);
					found = 1;
					break;
				}
			}
			if (!found) {
				printf("%.4lf\t", 0.0);
			}
		}
		printf("\n");
	}
}

CSRMatrix poisson1D_csr_matrix(int nb_equations) {
	CSRMatrix csr;
	csr.nb_rows = nb_equations;
	csr.nb_cols = nb_equations;
	csr.nb_non_zero = 3 * nb_equations - 2;
	csr.values = (double*)malloc(csr.nb_non_zero * sizeof(double));
	csr.col_index = (int*)malloc(csr.nb_non_zero * sizeof(int));
	csr.row_ptr = (int*)malloc((nb_equations + 1) * sizeof(int));
	// first line has 2 and -1
	csr.values[0] = 2;
	csr.col_index[0] = 0;
	csr.values[1] = -1;
	csr.col_index[1] = 1;
	csr.row_ptr[0] = 0;
	csr.row_ptr[1] = 2;
	// middle lines have -1, 2, -1
	for (int i = 1; i < nb_equations - 1; i++) {
		csr.values[3 * i - 1] = -1;
		csr.col_index[3 * i - 1] = i - 1;
		csr.values[3 * i] = 2;
		csr.col_index[3 * i] = i;
		csr.values[3 * i + 1] = -1;
		csr.col_index[3 * i + 1] = i + 1;
		csr.row_ptr[i + 1] = 3 * i + 2;
	}
	// last line has -1 and 2
	csr.values[3 * nb_equations - 3] = -1;
	csr.col_index[3 * nb_equations - 3] = nb_equations - 2;
	csr.values[3 * nb_equations - 2] = 2;
	csr.col_index[3 * nb_equations - 2] = nb_equations - 1;
	csr.row_ptr[nb_equations] = 3 * nb_equations;
	return csr;
}

double* csc_to_dense_col_major(CSCMatrix* csc) {
	double* dense = (double*)malloc(csc->nb_rows * csc->nb_cols * sizeof(double));
	memset(dense, 0, csc->nb_rows * csc->nb_cols * sizeof(double));
	for (int i = 0; i < csc->nb_cols; i++) {
		for (int j = csc->col_ptr[i]; j < csc->col_ptr[i + 1]; j++) {
			dense[csc->row_index[j] + i * csc->nb_rows] = csc->values[j];
		}
	}
	return dense;
}

CSCMatrix dense_col_major_to_csc(double* dense, int nb_rows, int nb_cols) {
	CSCMatrix csc;
	csc.nb_rows = nb_rows;
	csc.nb_cols = nb_cols;
	csc.nb_non_zero = 0;
	for (int i = 0; i < nb_cols; i++) {
		for (int j = 0; j < nb_rows; j++) {
			if (dense[j + i * nb_rows] != 0) {
				csc.nb_non_zero++;
			}
		}
	}
	csc.values = (double*)malloc(csc.nb_non_zero * sizeof(double));
	csc.row_index = (int*)malloc(csc.nb_non_zero * sizeof(int));
	csc.col_ptr = (int*)malloc((nb_cols + 1) * sizeof(int));
	int non_zero_idx = 0;
	for (int i = 0; i < nb_cols; i++) {
		csc.col_ptr[i] = non_zero_idx;
		for (int j = 0; j < nb_rows; j++) {
			if (dense[j + i * nb_rows] != 0) {
				csc.values[non_zero_idx] = dense[j + i * nb_rows];
				csc.row_index[non_zero_idx] = j;
				non_zero_idx++;
			}
		}
	}
	csc.col_ptr[nb_cols] = non_zero_idx;
	return csc;
}

void print_csc_matrix(CSCMatrix* csc) {
	for (int i = 0; i < csc->nb_rows; i++) {
		for (int j = 0; j < csc->nb_cols; j++) {
			int found = 0;
			for (int k = csc->col_ptr[j]; k < csc->col_ptr[j + 1]; k++) {
				if (csc->row_index[k] == i) {
					printf("%lf\t", csc->values[k]);
					found = 1;
					break;
				}
			}
			if (!found) {
				printf("%lf\t", 0.0);
			}
		}
		printf("\n");
	}
}

CSCMatrix poisson1D_csc_matrix(int nb_equations) {
	CSCMatrix csc;
	csc.nb_rows = nb_equations;
	csc.nb_cols = nb_equations;
	csc.nb_non_zero = 3 * nb_equations - 2;
	csc.values = (double*)malloc(csc.nb_non_zero * sizeof(double));
	csc.row_index = (int*)malloc(csc.nb_non_zero * sizeof(int));
	csc.col_ptr = (int*)malloc((nb_equations + 1) * sizeof(int));
	// first column has 2 and -1
	csc.values[0] = 2;
	csc.row_index[0] = 0;
	csc.values[1] = -1;
	csc.row_index[1] = 1;
	csc.col_ptr[0] = 0;
	csc.col_ptr[1] = 2;
	// middle columns have -1, 2, -1
	for (int i = 1; i < nb_equations - 1; i++) {
		csc.values[3 * i - 1] = -1;
		csc.row_index[3 * i - 1] = i - 1;
		csc.values[3 * i] = 2;
		csc.row_index[3 * i] = i;
		csc.values[3 * i + 1] = -1;
		csc.row_index[3 * i + 1] = i + 1;
		csc.col_ptr[i + 1] = 3 * i + 2;
	}
	// last column has -1 and 2
	csc.values[3 * nb_equations - 3] = -1;
	csc.row_index[3 * nb_equations - 3] = nb_equations - 2;
	csc.values[3 * nb_equations - 2] = 2;
	csc.row_index[3 * nb_equations - 2] = nb_equations - 1;
	csc.col_ptr[nb_equations] = 3 * nb_equations;
	return csc;
}

CSRMatrix csr_clone(CSRMatrix* csr) {
	CSRMatrix clone;
	clone.nb_rows = csr->nb_rows;
	clone.nb_cols = csr->nb_cols;
	clone.nb_non_zero = csr->nb_non_zero;
	clone.values = (double*)malloc(clone.nb_non_zero * sizeof(double));
	clone.col_index = (int*)malloc(clone.nb_non_zero * sizeof(int));
	clone.row_ptr = (int*)malloc((clone.nb_rows + 1) * sizeof(int));
	memcpy(clone.values, csr->values, clone.nb_non_zero * sizeof(double));
	memcpy(clone.col_index, csr->col_index, clone.nb_non_zero * sizeof(int));
	memcpy(clone.row_ptr, csr->row_ptr, (clone.nb_rows + 1) * sizeof(int));
	return clone;
}

CSCMatrix csc_clone(CSCMatrix* csc) {
	CSCMatrix clone;
	clone.nb_rows = csc->nb_rows;
	clone.nb_cols = csc->nb_cols;
	clone.nb_non_zero = csc->nb_non_zero;
	clone.values = (double*)malloc(clone.nb_non_zero * sizeof(double));
	clone.row_index = (int*)malloc(clone.nb_non_zero * sizeof(int));
	clone.col_ptr = (int*)malloc((clone.nb_cols + 1) * sizeof(int));
	memcpy(clone.values, csc->values, clone.nb_non_zero * sizeof(double));
	memcpy(clone.row_index, csc->row_index, clone.nb_non_zero * sizeof(int));
	memcpy(clone.col_ptr, csc->col_ptr, (clone.nb_cols + 1) * sizeof(int));
	return clone;
}

void csr_free(CSRMatrix* csr) {
	free(csr->values);
	free(csr->col_index);
	free(csr->row_ptr);
}

void csc_free(CSCMatrix* csc) {
	free(csc->values);
	free(csc->row_index);
	free(csc->col_ptr);
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

// TODO
void dcscmv(char is_M_transposed, double alpha, CSCMatrix* M, double* x, size_t incx, double beta, double* y, size_t incy) {
	// If M is not transposed, we can directly do the dot product of each row of M with x and not have to use extra memory
	if (is_M_transposed == 'N') {
		for (int i = 0; i < M->nb_cols; i++) {
			double sum = 0;
			for (int nz_idx = M->col_ptr[i]; nz_idx < M->col_ptr[i + 1]; nz_idx++) {
				sum += M->values[nz_idx] * x[M->row_index[nz_idx] * incx];
			}
			y[i * incy] = alpha * sum + beta * y[i * incy];
		}
	}
	// If M is transposed, we can't do it line by line, we need to store each result in their appropriate column and then sum them up
	else {
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
}

double csr_elem_at(CSRMatrix* csr, int i, int j) {
	for (int nz_idx = csr->row_ptr[i]; nz_idx < csr->row_ptr[i + 1]; nz_idx++) {
		if (csr->col_index[nz_idx] == j) {
			return csr->values[nz_idx];
		}
	}
	return 0;
}

double csc_elem_at(CSCMatrix* csc, int i, int j) {
	for (int nz_idx = csc->col_ptr[j]; nz_idx < csc->col_ptr[j + 1]; nz_idx++) {
		if (csc->row_index[nz_idx] == i) {
			return csc->values[nz_idx];
		}
	}
	return 0;
}

CSRMatrix csr_from_diag(double* diag, int nb_elements) {
	CSRMatrix csr;
	csr.nb_rows = nb_elements;
	csr.nb_cols = nb_elements;
	csr.nb_non_zero = nb_elements;
	csr.values = (double*)malloc(csr.nb_non_zero * sizeof(double));
	csr.col_index = (int*)malloc(csr.nb_non_zero * sizeof(int));
	csr.row_ptr = (int*)malloc((nb_elements + 1) * sizeof(int));
	for (int i = 0; i < nb_elements; i++) {
		csr.values[i] = diag[i];
		csr.col_index[i] = i;
		csr.row_ptr[i] = i;
	}
	csr.row_ptr[nb_elements] = nb_elements;
	return csr;
}

CSCMatrix csc_from_diag(double* diag, int nb_elements) {
	CSCMatrix csc;
	csc.nb_rows = nb_elements;
	csc.nb_cols = nb_elements;
	csc.nb_non_zero = nb_elements;
	csc.values = (double*)malloc(csc.nb_non_zero * sizeof(double));
	csc.row_index = (int*)malloc(csc.nb_non_zero * sizeof(int));
	csc.col_ptr = (int*)malloc((nb_elements + 1) * sizeof(int));
	for (int i = 0; i < nb_elements; i++) {
		csc.values[i] = diag[i];
		csc.row_index[i] = i;
		csc.col_ptr[i] = i;
	}
	csc.col_ptr[nb_elements] = nb_elements;
	return csc;
}

CSRMatrix csr_from_tridiag(double* diag, double* lower, double* upper, int nb_elements) {
	CSRMatrix csr;
	csr.nb_rows = nb_elements;
	csr.nb_cols = nb_elements;
	csr.nb_non_zero = 3 * nb_elements - 2;
	csr.values = (double*)malloc(csr.nb_non_zero * sizeof(double));
	csr.col_index = (int*)malloc(csr.nb_non_zero * sizeof(int));
	csr.row_ptr = (int*)malloc((nb_elements + 1) * sizeof(int));
	// first line has 2 and -1
	csr.values[0] = diag[0];
	csr.col_index[0] = 0;
	csr.values[1] = upper[0];
	csr.col_index[1] = 1;
	csr.row_ptr[0] = 0;
	csr.row_ptr[1] = 2;
	// middle lines have -1, 2, -1
	for (int i = 1; i < nb_elements - 1; i++) {
		csr.values[3 * i - 1] = lower[i - 1];
		csr.col_index[3 * i - 1] = i - 1;
		csr.values[3 * i] = diag[i];
		csr.col_index[3 * i] = i;
		csr.values[3 * i + 1] = upper[i];
		csr.col_index[3 * i + 1] = i + 1;
		csr.row_ptr[i + 1] = 3 * i + 2;
	}
	// last line has -1 and 2
	csr.values[3 * nb_elements - 3] = lower[nb_elements - 2];
	csr.col_index[3 * nb_elements - 3] = nb_elements - 2;
	csr.values[3 * nb_elements - 2] = diag[nb_elements - 1];
	csr.col_index[3 * nb_elements - 2] = nb_elements - 1;
	csr.row_ptr[nb_elements] = 3 * nb_elements;
	return csr;
}

CSCMatrix csc_from_tridiag(double* diag, double* lower, double* upper, int nb_elements) {
	CSCMatrix csc;
	csc.nb_rows = nb_elements;
	csc.nb_cols = nb_elements;
	csc.nb_non_zero = 3 * nb_elements - 2;
	csc.values = (double*)malloc(csc.nb_non_zero * sizeof(double));
	csc.row_index = (int*)malloc(csc.nb_non_zero * sizeof(int));
	csc.col_ptr = (int*)malloc((nb_elements + 1) * sizeof(int));
	// first column has 2 and -1
	csc.values[0] = diag[0];
	csc.row_index[0] = 0;
	csc.values[1] = upper[0];
	csc.row_index[1] = 1;
	csc.col_ptr[0] = 0;
	csc.col_ptr[1] = 2;
	// middle columns have -1, 2, -1
	for (int i = 1; i < nb_elements - 1; i++) {
		csc.values[3 * i - 1] = lower[i - 1];
		csc.row_index[3 * i - 1] = i - 1;
		csc.values[3 * i] = diag[i];
		csc.row_index[3 * i] = i;
		csc.values[3 * i + 1] = upper[i];
		csc.row_index[3 * i + 1] = i + 1;
		csc.col_ptr[i + 1] = 3 * i + 2;
	}
	// last column has -1 and 2
	csc.values[3 * nb_elements - 3] = lower[nb_elements - 2];
	csc.row_index[3 * nb_elements - 3] = nb_elements - 2;
	csc.values[3 * nb_elements - 2] = diag[nb_elements - 1];
	csc.row_index[3 * nb_elements - 2] = nb_elements - 1;
	csc.col_ptr[nb_elements] = 3 * nb_elements;
	return csc;
}

CSRMatrix csr_from_lower_triangular(double* mat, int nb_rows) {
	// mat is a lower triangular matrix, we can directly copy the values
	CSRMatrix csr;
	csr.nb_rows = nb_rows;
	csr.nb_cols = nb_rows;
	csr.nb_non_zero = nb_rows * (nb_rows + 1) / 2;
	csr.values = (double*)malloc(csr.nb_non_zero * sizeof(double));
	csr.col_index = (int*)malloc(csr.nb_non_zero * sizeof(int));
	csr.row_ptr = (int*)malloc((nb_rows + 1) * sizeof(int));
	int non_zero_idx = 0;
	int current_idx = 0;
	for (int i = 0; i < nb_rows; i++) {
		csr.row_ptr[i] = non_zero_idx;
		for (int j = 0; j <= i; j++) {
			csr.values[non_zero_idx] = mat[current_idx];
			csr.col_index[non_zero_idx] = j;
			non_zero_idx++;
			current_idx++;
		}
	}
	csr.row_ptr[nb_rows] = non_zero_idx;
	return csr;
}

CSCMatrix csc_from_lower_triangular(double* mat, int nb_rows) {
	// TODO
	exit(1);
}