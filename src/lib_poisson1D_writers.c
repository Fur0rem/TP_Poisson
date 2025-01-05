/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename) {
	FILE* file;
	int ii, jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (ii = 0; ii < (*lab); ii++) {
			for (jj = 0; jj < (*la); jj++) {
				fprintf(file, "%lf\t", AB[ii * (*la) + jj]);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
	else {
		perror(filename);
	}
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename) {
	FILE* file;
	int ii, jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (ii = 0; ii < (*la); ii++) {
			for (jj = 0; jj < (*lab); jj++) {
				fprintf(file, "%lf\t", AB[ii * (*lab) + jj]);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
	else {
		perror(filename);
	}
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename) {
	FILE* file;
	int jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 1; jj < (*la); jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj, jj + 1, AB[(*la) + jj]);
		}
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj + 1, jj + 1, AB[2 * (*la) + jj]);
		}
		for (jj = 0; jj < (*la) - 1; jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj + 2, jj + 1, AB[3 * (*la) + jj]);
		}
		fclose(file);
	}
	else {
		perror(filename);
	}
}

void write_vec(double* vec, int* la, char* filename) {
	int jj;
	FILE* file;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%lf\n", vec[jj]);
		}
		fclose(file);
	}
	else {
		perror(filename);
	}
}

// TODO: clean this up
void plot_convergence_history(double* vec, int size, char* filename) {
	// Use matplotlib to plot the convergence history
	char command[10000];
	char* filename_without_convergence = filename + strlen("rapport/convergence_");
	char* filename_without_convergence_clone = strdup(filename_without_convergence);

	// Put upper case to the first letter and every letter after a space, dash or underscore
	filename_without_convergence_clone[0] = filename_without_convergence_clone[0] - ('a' - 'A');
	for (int i = 1; i < strlen(filename_without_convergence_clone); i++) {
		if (filename_without_convergence_clone[i] == ' ' || filename_without_convergence_clone[i] == '-' ||
		    filename_without_convergence_clone[i] == '_') {
			if (filename_without_convergence_clone[i] == '_') {
				filename_without_convergence_clone[i] = ' ';
			}
			if (filename_without_convergence_clone[i + 1] >= 'a' && filename_without_convergence_clone[i + 1] <= 'z') {
				filename_without_convergence_clone[i + 1] = filename_without_convergence_clone[i + 1] - ('a' - 'A');
			}
		}
	}

	printf("Plotting the convergence history of the method %s\n", filename_without_convergence_clone);

	sprintf(command, "python3 -c \"import matplotlib.pyplot as plt; import numpy as np;");
	// Add the data
	char data[10000] = {0};
	sprintf(data, "data = np.array([");
	for (int i = 0; i < size; i++) {
		char temp[100] = {0};
		sprintf(temp, "[%d, %lf],", i, vec[i]);
		strcat(data, temp);
	}
	strcat(data, "]);");
	// Add the plot
	char plot[1000] = {0};
	sprintf(plot,
		"plt.plot(data[:,0], data[:,1]); plt.title('Convergence de la méthode de %s'); plt.xlabel('Nombre d\\'itérations'); "
		"plt.ylabel('Erreur avant relative'); plt.savefig('%s.svg')\"",
		filename_without_convergence_clone,
		filename);
	strcat(command, data);
	strcat(command, plot);
	int res = system(command);
	if (res) {
		perror("Error while plotting the convergence history");
	}
}

void write_xy(double* vec, double* x, int* la, char* filename) {
	int jj;
	FILE* file;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
		}
		fclose(file);
	}
	else {
		perror(filename);
	}
}
