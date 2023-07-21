#ifndef SSPD3D_H_INCLUDED
#define SSPD3D_H_INCLUDED

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
#include "moreMath.h"


float point_to_seg(float* p, float* s1, float* s2); //calcula la menor dista
float point_to_seg_2d(float* p, float* s1, float* s2);
float point_to_fiber(float* p, fiber f, float (*pnt2seg)(float*,float*,float*)); //Entrega la menor distancia entre un punto y cualquier otro punto perteneciente a una curva(fibra)
float e_spd(fiber f1, fiber f2,float (*pnt2seg)(float*,float*,float*));//Calcula la distancia "point_to_fiber" promedio entre todos los puntos d una fibra con respecto a otra
float e_sspd(fiber f1, fiber f2,float (*pnt2seg)(float*,float*,float*)); //Hace el calculo anterior, con ambas fibras y obtiene un promedio (DSSDPD)
void ssdpd_matrix(fasciculus b1,fasciculus b2,float matrix[b1.numFibers][b2.numFibers],float (*pnt2seg)(float*,float*,float*));// Calcula todas las distancias entre todas las fibras de 2 fasciculus
float** OMP_ssdpd_matrix_op(fasciculus b1,float (*pnt2seg)(float*,float*,float*));//Versi�n optimizada del codigo anterior cuando se usa con el mismo fasciculus

// Adaptaciones de Python
float point_to_seg_py(float* p, float* s1, float* s2, float dps1, float dps2, float dss);
float point_to_seg_2d_py(float* p, float* s1, float* s2, float dps1, float dps2, float dss);
float point_to_fiber_py(float* p, fiber f, float mdist_p[f.numPts], float t_dist[f.numPts-1], int dim);
float e_spd_py(fiber f1, fiber f2, float mdist[f1.numPts][f2.numPts], float t2_dist[f2.numPts-1], int dim);
float e_sspd_py(fiber f1, fiber f2, int dim);
float** ssdpd_matrix_py(fasciculus b1, fasciculus b2, int dim);
float** ssdpd_matrix_py_op(fasciculus b1, int dim);

#endif // SSPD3D_H_INCLUDED
