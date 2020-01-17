#include <conio.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <omp.h>

FILE *fl;
errno_t err;
const double eps = 0.00000001;

void fill_random(double **&arr, double *&vector, int n) {
	arr = new double*[n];
	for (int i = 0; i < n; i++)
		arr[i] = new double[n];
	srand(time(0));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				arr[i][j] = rand() % 10000 - rand() % 10000; //Заполнение массива случайными числами от -9999 до 9999
				//printf_s("%.2lf ", arr[i][j]);// - Проверял случайность чисел и, чтобы иметь данные необходимые для дальнейшей проверки умножений.
			}
			//printf_s("\n");
		}
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = rand() % 10000 - rand() % 10000;
	//	printf_s("%.2lf ", vector[i]);
	}
}

void write_time(double search_time, int n, const char *name_file) {
	err = fopen_s(&fl, *&name_file, "a");
	if (err == 0) {
		printf_s("%d \n", n);
	}
	else {
		printf_s("Ошибка открытия");
		return;
	}
	fprintf_s(fl, "Время работы: %f - %d элементов\n", search_time, n);
	fclose(fl);
}

double *Gaus(int n, double **arr, double *vector) {
	int k = 0, g = 0;
	double start = omp_get_wtime();
	for (k = 0; k < n; k++) {
		for (int i = k; i < n; i++) {
			g = 0;
			double t = arr[i][k];
			if (abs(t) < eps) continue; // Нулевой коэффициент, следовательно делить и вычитать не нужно.
			for (int j = 0; j < n; j++) {  //Делим на диагональный элемент соответствующей строки.
				arr[i][j] = arr[i][j] / t;
			}
			vector[i] = vector[i] / t;
			if (i == k) continue; // Чтобы не вычитать уравнения из себя же.
			for (int j = 0; j < n; j++) { // Вычитаю k-ую строку из i-ой.
				arr[i][j] -= arr[k][j];
				if (abs(arr[i][j]) < eps) {
					arr[i][j] = 0;
				}
			}
			vector[i] -= vector[k];
		}
	}
	double *x = NULL;
	x = new double[n];
	int free_var = 1; // Для частного решения при наличии свободных переменных.
	for (k = n - 1; k >= 0; k--) {
		g = 0;
		for (int j = k; j < n; j++) {
			if (abs(arr[k][j]) < eps) g++; // Проверка на нулевую строку (без вектора правой части).
		}
		if ((vector[k] == 0) && (g == n - k)) {
			x[k] = free_var;
			free_var++;
		}
		else if ((vector[k] != 0) && (g == n - k)) {
			printf_s("Система не имеет решений! ");
			write_time(omp_get_wtime() - start, n, "no_parall.txt");
			return 0;
		}
		else {
			x[k] = vector[k];
		}
		for (int i = 0; i < k; i++) {
			vector[i] = vector[i] - arr[i][k] * x[k];
		}
	}
	write_time(omp_get_wtime() - start, n, "no_parall.txt");
	return x;
}

double *Parall_Gaus(int n, double **arr, double *vector) {
	int k = 0, g = 0;
	double start = omp_get_wtime();
	for (k = 0; k < n; k++) {
	#pragma omp parallel for
		for (int i = k; i < n; i++) {
			g = 0;
			double t = arr[i][k];
			if (abs(t) < eps) continue; // Нулевой коэффициент, следовательно делить и вычитать не нужно.
			for (int j = 0; j < n; j++) {  //Делим на диагональный элемент соответствующей строки.
				arr[i][j] = arr[i][j] / t;
			}
			vector[i] = vector[i] / t;
			if (i == k) continue; // Чтобы не вычитать уравнения из себя же.
			for (int j = 0; j < n; j++) { // Вычитаю k-ую строку из i-ой.
				arr[i][j] -= arr[k][j];
				if (abs(arr[i][j]) < eps) {
					arr[i][j] = 0;
				}
			}
			vector[i] -= vector[k];
		}
	}
	double *x = NULL;
	x = new double[n];
	int free_var = 1; // Для частного решения при наличии свободных переменных.
	for (k = n - 1; k >= 0; k--) {
		g = 0;
		#pragma omp parallel for
		for (int j = k; j < n; j++) {
			if (abs(arr[k][j]) < eps) g++; // Проверка на нулевую строку (без вектора правой части).
		}
		if ((vector[k] == 0) && (g == n - k)) {
			x[k] = free_var;
			free_var++;
		}
		else if ((vector[k] != 0) && (g == n -k)) {
			printf_s("Система не имеет решений! ");
			write_time(omp_get_wtime()-start, n, "parall.txt");
			return 0;
		}
		else {
			x[k] = vector[k];
		}
		#pragma omp parallel for
		for (int i = 0; i < k; i++) {
			vector[i] = vector[i] - arr[i][k] * x[k];
		}
	}
	write_time(omp_get_wtime()-start, n, "parall.txt");
	return x;
}

void main() {
	setlocale(LC_CTYPE, "Rus");
	int n, k = 0;
	
	double **arr = NULL, *vector = NULL, *x = NULL;
	for (n = 300; n <= 1300; n += 100) {
		fill_random(arr, vector, n);
		//x = Gaus(n, arr, vector);
		x = Parall_Gaus(n, arr, vector);
		for (int i = 0; i < n; i++) {
			delete arr[i];
		}
		delete arr;
		delete[] vector;
	}
	printf_s("Job's done! ");
	_getch();

}