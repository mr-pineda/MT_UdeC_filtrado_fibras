#ifndef ISOMAP_H_INCLUDED
#define ISOMAP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

#define inf 1000000.0

#include "moreMath.h"

void dist_mat_op(float **coord, int num_points, float matrix[num_points][num_points]);
void swap_f(float *xp, float *yp);
void swap_i(int *xp, int *yp);

void find_kmin_index(int n, float arr[n], int k, int k_min[k]);
void knn_list(int n, float d_mat[n][n], int k, int knn_mat[n][k]);
bool n_in_array(int n, int dim, int arr[dim]);
void graph_mat(int n, float d_mat[n][n], int k, float graph[n][n]);

int find_min_index(int n, float arr[n]);
bool all_visited(int n, bool visited[n]);
bool dijkstra(int n, float graph[n][n], int start, float short_path[n]);
bool geodesic_mat(int n, float graph[n][n], float geo_dist[n][n]);
void fiber2geomat(int n, float **coord, int k, float geomat[n][n]);

void centering_matrix(int n, float c_mat[n][n]);
void sqr_dot_product(int dim, float M[dim][dim], float N[dim][dim], float X[dim][dim]);
void sqr_matrix(int dim, float M[dim][dim]);
void sqr_matrix_op(int dim, float M[dim][dim]);
void dbl_centered_matrix(int n, float mat[n][n], float B[n][n]);
void dbl_centered_matrix_op(int n, float mat[n][n], float B[n][n]);
void MDS_2D(int n, float **coord, int iter, float MDS[2][n], void (*dc_mat)(int, float[n][n], float[n][n]));
void ISOMAP_2D(int n, float **coord, int k, int iter, float isomap[2][n], void (*dc_mat)(int, float[n][n], float[n][n]));

float rad2deg(float rad);
float deg2rad(float deg);
void rotate_points(float rad, int n, float arr_2d[2][n]);
float rad_m(float p1[2], float p2[2]);
void flip_ud(int n, float arr_2d[2][n]);
void flip_lr(int n, float arr_2d[2][n]);

#endif //
