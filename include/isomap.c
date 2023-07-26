#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <omp.h>

#include "isomap.h"
#include "moreMath.h"
#include "jacobi_eigenvalue.h"

/// Algoritmos para obtener distancias y func b�sicas auxiliares

void dist_mat_op(float **coord, int num_points, float matrix[num_points][num_points])
{
    for (int i = 0; i < num_points; i++)
    {
        matrix[i][i] = 0;
        for (int j = i; j < num_points; j++)
        {
            matrix[i][j] = eucDisPP(coord[i], coord[j]);
            matrix[j][i] = matrix[i][j];
        }
    }
    return;
}

void swap_f(float *xp, float *yp)
{
    float temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap_i(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/// Obtenci�n de KNN y contrucci�n del grafo
void find_kmin_index(int n, float arr[n], int k, int k_min[k])
{
    // Adaptaci�n de l algoritmo de selection sort, encuentra los k valores minimos de un vector
    // y retorna sus �ndices.
    int i, j, min_idx;
    int idx_vec[n];
    float temp[n];

    for (i = 0; i < n; i++)
    {
        temp[i] = arr[i];
        idx_vec[i] = i;
    }

    for (i = 0; i < k; i++)
    {
        min_idx = i;
        for (j = i + 1; j < n; j++)
            if (temp[j] < temp[min_idx])
                min_idx = j;
        swap_f(&temp[min_idx], &temp[i]);
        swap_i(&idx_vec[min_idx], &idx_vec[i]);
        k_min[i] = idx_vec[i];
    }
    return;
}

void knn_list(int n, float d_mat[n][n], int k, int knn_mat[n][k])
{
    for (int i = 0; i < n; i++)
    {
        find_kmin_index(n, d_mat[i], k, knn_mat[i]);
    }
    return;
}

bool n_in_array(int n, int dim, int arr[dim])
{
    int i = 0;
    while (i < dim)
    {
        if (n == arr[i])
        {
            return true;
        }
        i++;
    }
    return false;
}

void graph_mat(int n, float d_mat[n][n], int k, float graph[n][n])
{
    int k_mat[n][k];
    knn_list(n, d_mat, k, k_mat);

    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
            graph[row][col] = inf;
    }
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (n_in_array(col, k, k_mat[row]))
            {
                graph[row][col] = d_mat[row][col];
                graph[col][row] = d_mat[row][col];
            }
        }
    }
    return;
}

/// Obtenci�n de ruta m�s corta (Algoritmo de dijkstra)

int find_min_index(int n, float arr[n])
{
    float min = arr[0];
    int min_idx = 0;
    for (int i = 1; i < n; i++)
    {
        if (min > arr[i])
        {
            min = arr[i];
            min_idx = i;
        }
    }
    return min_idx;
}

bool all_visited(int n, bool visited[n])
{
    for (int i = 0; i < n; i++)
    {
        if (!visited[i])
            return false;
    }
    return true;
}

bool all_infinites(int n, float cand[n])
{
    for (int i = 0; i < n; i++)
    {
        if (cand[i] != inf)
            return false;
    }
    return true;
}

bool dijkstra(int n, float graph[n][n], int start, float short_path[n])
{
    for (int i = 0; i < n; i++)
    {
        short_path[i] = inf;
    }
    bool visited[n];
    for (int i = 0; i < n; i++)
    {
        visited[i] = false;
    }
    short_path[start] = 0;

    while (!all_visited(n, visited))
    {
        float candidates[n];
        for (int i = 0; i < n; i++)
        {
            if (!visited[i])
                candidates[i] = short_path[i];
            else
                candidates[i] = inf;
        }

        if (all_infinites(n, candidates))
            return true;
        int m_idx = find_min_index(n, candidates);

        for (int i = 0; i < n; i++)
        {
            float new_dist = candidates[m_idx] + graph[m_idx][i];
            if (new_dist < short_path[i])
                short_path[i] = new_dist;
        }
        visited[m_idx] = true;
    }
    return false;
}

bool geodesic_mat(int n, float graph[n][n], float geo_dist[n][n])
{
    for (int i = 0; i < n; i++)
    {
        if (dijkstra(n, graph, i, geo_dist[i]))
            return true;
    }
    return false;
}

void fiber2geomat(int n, float **coord, int k, float geomat[n][n])
{
    float dist[n][n];
    float graph[n][n];
    int knn = k + 1;
    bool error = false;

    dist_mat_op(coord, n, dist);
    do
    {
        graph_mat(n, dist, knn, graph);
        error = geodesic_mat(n, graph, geomat);

        knn++;
    } while (error);

    return;
}

/// Aplicaci�n del algoritmo de MDS

void centering_matrix(int n, float c_mat[n][n])
{
    float trian = -1.0 / n;
    float diag = (n - 1.0) / n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                c_mat[i][j] = diag;
            else
                c_mat[i][j] = trian;
        }
    }

    return;
}

void sqr_dot_product(int dim, float M[dim][dim], float N[dim][dim], float X[dim][dim])
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            X[i][j] = 0;
            for (int k = 0; k < dim; k++)
            {
                X[i][j] += M[i][k] * N[k][j];
            }
        }
    }
    return;
}

void sqr_matrix(int dim, float M[dim][dim])
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            M[i][j] *= M[i][j];
        }
    }
}

void sqr_matrix_op(int dim, float M[dim][dim])
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = i + 1; j < dim; j++)
        {
            M[i][j] *= M[i][j];
            M[j][i] = M[i][j];
        }
    }
}

void dbl_centered_matrix(int n, float mat[n][n], float B[n][n])
{
    sqr_matrix(n, mat);

    float C[n][n];
    centering_matrix(n, C);

    float C2[n][n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            C2[i][j] = (-0.5) * C[i][j];
    }

    float mul1[n][n];
    sqr_dot_product(n, C2, mat, mul1);
    sqr_dot_product(n, mul1, C, B);

    return;
}

void dbl_centered_matrix_op(int n, float mat[n][n], float B[n][n])
{
    sqr_matrix_op(n, mat);

    float mean_col[n];
    float mean_row[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B[i][j] = mat[i][j];
        }
    }

    for (int i = 0; i < n; i++)
    {
        mean_col[i] = 0;
        for (int j = 0; j < n; j++)
        {
            mean_col[i] += B[j][i];
        }
        mean_col[i] /= n;
        for (int j = 0; j < n; j++)
        {
            B[j][i] -= mean_col[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        mean_row[i] = 0;
        for (int j = 0; j < n; j++)
        {
            mean_row[i] += B[i][j];
        }
        mean_row[i] /= n;
        for (int j = 0; j < n; j++)
        {
            B[i][j] -= mean_row[i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            B[i][j] *= (-0.5);
    }
    return;
}

void MDS_2D(int n, float **coord, int iter, float MDS[2][n], void (*dc_mat)(int, float[n][n], float[n][n]))
{
    float distmat[n][n];
    dist_mat_op(coord, n, distmat);

    float dcMat[n][n];
    dc_mat(n, distmat, dcMat);
    float eig_vec[n][n];
    float eig_val[n];

    int iter_num;
    int rot_num;

    jacobi_eigenvalue(n, dcMat, iter, eig_vec, eig_val, &iter_num, &rot_num);
    for (int i = 0; i < n; i++)
    {
        MDS[0][i] = eig_vec[0][i] * sqrt(eig_val[0]);
        MDS[1][i] = eig_vec[1][i] * sqrt(eig_val[1]);
    }
    return;
}

void ISOMAP_2D(int n, float **coord, int k, int iter, float isomap[2][n], void (*dc_mat)(int, float[n][n], float[n][n]))
{
    float geomat[n][n];
    fiber2geomat(n, coord, k, geomat);

    float dcMat[n][n];
    dc_mat(n, geomat, dcMat);

    float eig_vec[n][n];
    float eig_val[n];
    int iter_num;
    int rot_num;

    jacobi_eigenvalue(n, dcMat, iter, eig_vec, eig_val, &iter_num, &rot_num);

    for (int i = 0; i < n; i++)
    {
        isomap[0][i] = eig_vec[0][i] * sqrt(eig_val[0]);
        isomap[1][i] = eig_vec[1][i] * sqrt(eig_val[1]);
    }
    return;
}

float rad2deg(float rad)
{
    return rad * 180.0 / M_PI;
}

float deg2rad(float deg)
{
    return deg * M_PI / 180.0;
}

void rotate_points(float rad, int n, float arr_2d[2][n])
{
    float x, y;
    for (int i = 0; i < n; i++)
    {
        x = arr_2d[0][i] * cos(rad) - arr_2d[1][i] * sin(rad);
        y = arr_2d[0][i] * sin(rad) + arr_2d[1][i] * cos(rad);
        arr_2d[0][i] = x;
        arr_2d[1][i] = y;
    }
}

float rad_m(float p1[2], float p2[2])
{
    float m = (p2[1] - p1[1]) / (p2[0] - p1[0]);

    return atan(m);
};

void flip_ud(int n, float arr_2d[2][n])
{
    for (int i = 0; i < n; i++)
    {
        arr_2d[1][i] = -arr_2d[1][i];
    }
}

void flip_lr(int n, float arr_2d[2][n])
{
    for (int i = 0; i < n; i++)
    {
        arr_2d[0][i] = -arr_2d[0][i];
    }
}
