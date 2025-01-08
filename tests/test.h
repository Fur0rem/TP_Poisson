/**
 * @file test.h
 * @brief Test utilities
 */

#ifndef TEST_H
#define TEST_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

/**
 * @brief Print a matrix in column major format
 */
void print_mat_col_major(double* mat, int rows, int cols) {
        for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                        int idx = j * rows + i;
                        printf("%d : %lf\t", idx, mat[j * rows + i]);
                }
                printf("\n");
        }
}

/**
 * @brief Print a matrix in row major format
 */
void print_mat_row_major(double* mat, int rows, int cols) {
        for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                        printf("%lf\t", mat[j * rows + i]);
                }
                printf("\n");
        }
}

/**
 * @brief Check if two matrices are approximately equal
 */
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

/**
 * @brief Check if two vectors are approximately equal
 */
bool vec_equals(double* vec1, double* vec2, int size) {
        for (int i = 0; i < size; i++) {
                if (fabs(vec1[i] - vec2[i]) > 1e-6) {
                        return false;
                }
        }
        return true;
}

/**
 * @brief Flattens a 2D matrix into a 1D matrix row major format
 */
double* flattened_row_major(int rows, int cols, double mat[rows][cols]) {
        double* res = (double*)malloc(rows * cols * sizeof(double));
        for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                        res[j * rows + i] = mat[i][j];
                }
        }
        return res;
}

/**
 * @brief Flattens a 2D matrix into a 1D matrix column major format
 */
double* row_to_col_major(double* mat, int rows, int cols) {
        double* res = (double*)malloc(rows * cols * sizeof(double));
        for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                        res[j * rows + i] = mat[i * cols + j];
                }
        }
        return res;
}

#define RED   "\033[0;31m" ///< Red color
#define GREEN "\033[0;32m" ///< Green color
#define RESET "\033[0m"    ///< Reset color

/**
 * @brief Assert a condition
 */
#define ASSERT(condition)                                                                                                                  \
        if (!(condition)) {                                                                                                                \
                printf(RED "Assertion failed: " #condition "\n" RESET "File: %s, Line: %d\n", __FILE__, __LINE__);                         \
                exit(1);                                                                                                                   \
        }                                                                                                                                  \
        else {                                                                                                                             \
                printf(GREEN "Assertion passed: " #condition "\n" RESET);                                                                  \
        }

#endif // TEST_H