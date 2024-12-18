#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
	double* average_nanos;
	size_t size;
	size_t capacity;
} BenchmarkResultList;

BenchmarkResultList* benchmark_result_list_new() {
	BenchmarkResultList* list = (BenchmarkResultList*)malloc(sizeof(BenchmarkResultList));
	list->average_nanos = (double*)malloc(sizeof(double) * 10);
	list->size = 0;
	list->capacity = 10;
	return list;
}

void benchmark_result_list_push(BenchmarkResultList* list, double average_nanos) {
	if (list->size == list->capacity) {
		list->capacity *= 2;
		list->average_nanos = (double*)realloc(list->average_nanos, sizeof(double) * list->capacity);
	}
	list->average_nanos[list->size] = average_nanos;
	list->size++;
}

void benchmark_result_list_free(BenchmarkResultList* list) {
	free(list->average_nanos);
	free(list);
}

void write_bench_results(BenchmarkResultList* list, char* filename) {
	char filename_with_extension_dat[1000];
	char filename_with_extension_svg[1000];
	char filename_without_benchmark[1000]; // Remove the benchmark_ prefix
	sprintf(filename_without_benchmark, "%s", filename + 10);

	sprintf(filename_with_extension_dat, "rapport/%s.dat", filename);
	sprintf(filename_with_extension_svg, "rapport/%s.svg", filename);

	FILE* file = fopen(filename_with_extension_dat, "w");
	for (int i = 0; i < list->size; i++) {
		size_t nbpoints = 100 + i * 100;
		fprintf(file, "%zu %f\n", nbpoints, list->average_nanos[i]);
	}
	fclose(file);

	// Plot the file using matplotlib and save it as svg
	char command[10000];
	sprintf(command,
		"python3 -c\
		\"import matplotlib.pyplot as plt;\
		import numpy as np;\
		data = np.loadtxt('%s');\
		plt.plot(data[:,0], data[:,1]);\
		plt.xlabel('Nombre de points');\
		plt.ylabel('Temps moyen (nanosecondes)');\
		plt.title('Résultat benchmark %s');\
		plt.savefig('%s')\"",
		filename_with_extension_dat,
		filename_without_benchmark,
		filename_with_extension_svg);
	system(command);
}

void benchmark() {
	BenchmarkResultList* list_trf = benchmark_result_list_new();
	BenchmarkResultList* list_tri = benchmark_result_list_new();
	BenchmarkResultList* list_sv = benchmark_result_list_new();

	double T0 = -5.0;
	double T1 = 5.0;

	for (int nbpoints = 100; nbpoints < 100000; nbpoints += 100) {
		printf("nbpoints: %d\n", nbpoints);
		int la = nbpoints - 2;
		int kv = 1;
		int ku = 1;
		int kl = 1;
		int lab = kv + kl + ku + 1;

		int info = 1;
		int NRHS = 1;

		double* AB = (double*)malloc(sizeof(double) * lab * la);
		double* RHS = (double*)malloc(sizeof(double) * la);
		double* EX_SOL = (double*)malloc(sizeof(double) * la);
		double* X = (double*)malloc(sizeof(double) * la);
		int* ipiv = (int*)calloc(la, sizeof(int));

		set_grid_points_1D(X, &la);
		set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
		set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

		const int NB_SAMPLES = 20;
		clock_t total_trf_and_trs = 0;
		clock_t total_sv = 0;
		for (int i = 0; i < NB_SAMPLES; i++) {
			set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
			clock_t start = clock();
			dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
			dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
			clock_t end = clock();
			total_trf_and_trs += end - start;

			set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
			start = clock();
			dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
			end = clock();
			total_sv += end - start;
		}

		double average_trf_and_trs = (double)total_trf_and_trs / (double)NB_SAMPLES;
		double average_sv = (double)total_sv / (double)NB_SAMPLES;
		benchmark_result_list_push(list_trf, average_trf_and_trs);
		benchmark_result_list_push(list_sv, average_sv);

		free(AB);
		free(RHS);
		free(EX_SOL);
		free(X);
		free(ipiv);
	}

	printf("Writing results\n");

	write_bench_results(list_trf, "benchmark_dgbtrf+dgbtrs");
	write_bench_results(list_sv, "benchmark_dgbsv");

	benchmark_result_list_free(list_trf);
	benchmark_result_list_free(list_sv);
}

int main() {
	benchmark();
	return 0;
}