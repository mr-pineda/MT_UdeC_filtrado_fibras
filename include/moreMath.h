#ifndef MOREMATH_H_INCLUDED
#define MOREMATH_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>

float eucDisPP(float *p1, float *p2);                              // Calcula la distancia euclidea entre 2 puntos en un espacio 3d
float eucDisPP_2d(float *p1, float *p2);                           // Calcula la distancia euclidea entre 2 puntos en un espacio 3d
float *matrix_sum(float **matrix, int rows, int cols, short axis); // Suma las filas o columnas de una matriz
float *matrix_sum_op(float **matrix, int len);                     // versi�n optimizada del codigo anterior
void m_transpose(int row, int col, float matrix[row][col], float mat_tr[col][row]);
void matrix_sum_py(int len, float **matrix, float sum[len]);
void quicksort(float *arr, int N);                      // algoritmo quicksor para ordenar valores
int *percentile(float arr[], int arr_size, float perc); // se obtiene los indices de los elementos que cumplen cierto percentil.
void freeMatrix(float **matrix, int rows);

#endif // MOREMATH_H_INCLUDED
