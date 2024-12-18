#ifndef TEST_H
#define TEST_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

void print_mat_col_major(double* mat, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int idx = j * rows + i;
			printf("%d : %lf\t", idx, mat[j * rows + i]);
		}
		printf("\n");
	}
}

void print_mat_row_major(double* mat, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%lf\t", mat[j * rows + i]);
		}
		printf("\n");
	}
}

bool mat_equals(double* mat1, double* mat2, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (fabs(mat1[j * rows + i] - mat2[j * rows + i]) > 1e-6) {
				return false;
			}
		}
	}
	return true;
}

double* flattened_row_major(int rows, int cols, double mat[rows][cols]) {
	double* res = (double*)malloc(rows * cols * sizeof(double));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			res[j * rows + i] = mat[i][j];
		}
	}
	return res;
}

double* row_to_col_major(double* mat, int rows, int cols) {
	double* res = (double*)malloc(rows * cols * sizeof(double));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			res[j * rows + i] = mat[i * cols + j];
		}
	}
	return res;
}

#define RED   "\033[0;31m"
#define GREEN "\033[0;32m"
#define RESET "\033[0m"

#define ASSERT(condition)                                                                                                                  \
	if (!(condition)) {                                                                                                                \
		printf(RED "Assertion failed: " #condition "\n" RESET "File: %s, Line: %d\n", __FILE__, __LINE__);                         \
		exit(1);                                                                                                                   \
	}                                                                                                                                  \
	else {                                                                                                                             \
		printf(GREEN "Assertion passed: " #condition "\n" RESET);                                                                  \
	}

#endif // TEST_H