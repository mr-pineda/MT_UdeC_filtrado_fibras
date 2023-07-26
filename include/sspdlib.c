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

#include "fiberlib.h"
#include "sspdlib.h"
#include "moreMath.h"

float point_to_seg(float *p, float *s1, float *s2)
{
    if (s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2])
        return eucDisPP(p, s1);

    float u, i[3], diff[3];

    diff[0] = s2[0] - s1[0];
    diff[1] = s2[1] - s1[1];
    diff[2] = s2[2] - s1[2];
    u = (p[0] - s1[0]) * diff[0] + (p[1] - s1[1]) * diff[1] + (p[2] - s1[2]) * diff[2];
    u = u / (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    if (u <= 0)
    {
        return eucDisPP(p, s1);
    }
    else if (u >= 1)
    {
        return eucDisPP(p, s2);
    }
    else
    {
        i[0] = s1[0] + u * diff[0];
        i[1] = s1[1] + u * diff[1];
        i[2] = s1[2] + u * diff[2];
        return eucDisPP(p, i);
    }
}

float point_to_seg_2d(float *p, float *s1, float *s2)
{
    if (s1[0] == s2[0] && s1[1] == s2[1])
        return eucDisPP_2d(p, s1);

    float u, i[2], diff[2];

    diff[0] = s2[0] - s1[0];
    diff[1] = s2[1] - s1[1];
    u = (p[0] - s1[0]) * diff[0] + (p[1] - s1[1]) * diff[1];
    u = u / (diff[0] * diff[0] + diff[1] * diff[1]);
    if (u <= 0)
    {
        return eucDisPP_2d(p, s1);
    }
    else if (u >= 1)
    {
        return eucDisPP_2d(p, s2);
    }
    else
    {
        i[0] = s1[0] + u * diff[0];
        i[1] = s1[1] + u * diff[1];
        return eucDisPP_2d(p, i);
    }
}

float point_to_fiber(float *p, fiber f, float (*pnt2seg)(float *, float *, float *))
{
    float dpf = pnt2seg(p, f.cordinates[0], f.cordinates[1]);
    float temp;
    for (int j = 2; j < f.numPts; j++)
    {
        temp = pnt2seg(p, f.cordinates[j - 1], f.cordinates[j]);
        dpf = (dpf < temp) ? dpf : temp;
    }
    return dpf;
}

float e_spd(fiber f1, fiber f2, float (*pnt2seg)(float *, float *, float *))
{
    float spd = 0;
    for (int i = 0; i < f1.numPts; i++)
    {
        spd += point_to_fiber(f1.cordinates[i], f2, pnt2seg);
    }
    return spd / f1.numPts;
}

float e_sspd(fiber f1, fiber f2, float (*pnt2seg)(float *, float *, float *))
{
    return (e_spd(f1, f2, pnt2seg) + e_spd(f2, f1, pnt2seg)) / 2;
}

void ssdpd_matrix(fasciculus b1, fasciculus b2, float matrix[b1.numFibers][b2.numFibers], float (*pnt2seg)(float *, float *, float *))
{
    for (int i = 0; i < b1.numFibers; i++)
    {
        for (int j = 0; j < b2.numFibers; j++)
        {
            matrix[i][j] = e_sspd(b1.fibers[i], b2.fibers[j], pnt2seg);
        }
    }
    return;
}

float **OMP_ssdpd_matrix_op(fasciculus b1, float (*pnt2seg)(float *, float *, float *))
{
    float **matrix = (float **)malloc(sizeof(float *) * b1.numFibers);
    for (int i = 0; i < b1.numFibers; i++)
    {
        matrix[i] = (float *)malloc(sizeof(float) * b1.numFibers);
        matrix[i][i] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < b1.numFibers; i++)
    {
        for (int j = i + 1; j < b1.numFibers; j++)
        {
            matrix[i][j] = e_sspd(b1.fibers[i], b1.fibers[j], pnt2seg);
            matrix[j][i] = matrix[i][j];
        }
    }
    return matrix;
}

/// Adaptaciones de Python usando memoria est�tica

float point_to_seg_py(float *p, float *s1, float *s2, float dps1, float dps2, float dss)
{

    if (s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2])
        return dps1;

    float u, i[3], diff[3];

    diff[0] = s2[0] - s1[0];
    diff[1] = s2[1] - s1[1];
    diff[2] = s2[2] - s1[2];
    u = (p[0] - s1[0]) * diff[0] + (p[1] - s1[1]) * diff[1] + (p[2] - s1[2]) * diff[2];
    // u = u/(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]); // m�s �ptimo
    u = u / (dss * dss);

    if (u <= 0)
    {
        return dps1;
    }
    else if (u >= 1)
    {
        return dps2;
    }
    else
    {
        i[0] = s1[0] + u * diff[0];
        i[1] = s1[1] + u * diff[1];
        i[2] = s1[2] + u * diff[2];
        return eucDisPP(p, i);
    }
}

float point_to_seg_2d_py(float *p, float *s1, float *s2, float dps1, float dps2, float dss)
{

    if (s1[0] == s2[0] && s1[1] == s2[1])
        return dps1;

    float u, i[2], diff[2];

    diff[0] = s2[0] - s1[0];
    diff[1] = s2[1] - s1[1];
    u = (p[0] - s1[0]) * diff[0] + (p[1] - s1[1]) * diff[1];
    // u = u/(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]); // m�s �ptimo
    u = u / (dss * dss);

    if (u <= 0)
    {
        return dps1;
    }
    else if (u >= 1)
    {
        return dps2;
    }
    else
    {
        i[0] = s1[0] + u * diff[0];
        i[1] = s1[1] + u * diff[1];
        return eucDisPP(p, i);
    }
}

float point_to_fiber_py(float *p, fiber f, float mdist_p[f.numPts], float t_dist[f.numPts - 1], int dim)
{
    float (*pnt2seg)(float *, float *, float *, float, float, float);
    switch (dim)
    {
    case 2:
        pnt2seg = &point_to_seg_2d_py;
        break;
    case 3:
        pnt2seg = &point_to_seg_py;
        break;
    default:
        printf("\n ERROR: Dimension %d no es valida.\n", dim);
        exit(-1);
    }

    float dpf = pnt2seg(p, f.cordinates[0], f.cordinates[1], mdist_p[0], mdist_p[1], t_dist[0]);
    float temp;
    for (int j = 2; j < f.numPts; j++)
    {
        temp = pnt2seg(p, f.cordinates[j - 1], f.cordinates[j], mdist_p[j - 1], mdist_p[j], t_dist[j - 1]);
        dpf = (dpf < temp) ? dpf : temp;
    }
    return dpf;
}

float e_spd_py(fiber f1, fiber f2, float mdist[f1.numPts][f2.numPts], float t2_dist[f2.numPts - 1], int dim)
{
    float spd = 0;
    for (int i = 0; i < f1.numPts; i++)
    {
        spd += point_to_fiber_py(f1.cordinates[i], f2, mdist[i], t2_dist, dim);
    }
    return spd / f1.numPts;
}

float e_sspd_py(fiber f1, fiber f2, int dim)
{

    float (*eudist)(float *, float *);

    switch (dim)
    {
    case 2:
        eudist = &eucDisPP_2d;
        break;
    case 3:
        eudist = &eucDisPP;
        break;
    default:
        printf("\n ERROR: Dimension %d no es valida.\n", dim);
        exit(-1);
    }

    float dist_mat[f1.numPts][f2.numPts];
    for (int i = 0; i < f1.numPts; i++)
    {
        dist_mat[i][i] = 0;
        for (int j = i + 1; j < f2.numPts; j++)
        {
            dist_mat[i][j] = eudist(f1.cordinates[i], f2.cordinates[j]);
            dist_mat[j][i] = dist_mat[i][j];
        }
    }
    float dist_mat_T[f2.numPts][f1.numPts];
    m_transpose(f1.numPts, f2.numPts, dist_mat, dist_mat_T);

    int t1_len = f1.numPts - 1;
    int t2_len = f2.numPts - 1;
    float t1_dist[t1_len], t2_dist[t2_len];
    for (int i = 0; i < t1_len; i++)
        t1_dist[i] = eudist(f1.cordinates[i], f1.cordinates[i + 1]);

    for (int i = 0; i < t2_len; i++)
        t2_dist[i] = eudist(f2.cordinates[i], f2.cordinates[i + 1]);

    return (e_spd_py(f1, f2, dist_mat, t2_dist, dim) + e_spd_py(f2, f1, dist_mat_T, t1_dist, dim)) / 2;
}

float **ssdpd_matrix_py(fasciculus b1, fasciculus b2, int dim)
{

    float **matrix = (float **)malloc(sizeof(float *) * b1.numFibers);
    for (int i = 0; i < b1.numFibers; i++)
    {
        matrix[i] = (float *)malloc(sizeof(float) * b1.numFibers);
    }

#pragma omp parallel for
    for (int i = 0; i < b1.numFibers; i++)
    {
        for (int j = 0; j < b1.numFibers; j++)
        {
            matrix[i][j] = e_sspd_py(b1.fibers[i], b2.fibers[j], dim);
        }
    }
    return matrix;
}

float **ssdpd_matrix_py_op(fasciculus b1, int dim)
{
    float **matrix = (float **)malloc(sizeof(float *) * b1.numFibers);
    for (int i = 0; i < b1.numFibers; i++)
    {
        matrix[i] = (float *)malloc(sizeof(float) * b1.numFibers);
        matrix[i][i] = 0;
    }
#pragma omp parallel
    {
#pragma omp for nowait
        for (int i = 0; i < b1.numFibers; i++)
        {
            for (int j = i + 1; j < b1.numFibers; j++)
            {
                matrix[i][j] = e_sspd_py(b1.fibers[i], b1.fibers[j], dim);
                matrix[j][i] = matrix[i][j];
            }
        }
    }

    return matrix;
}
