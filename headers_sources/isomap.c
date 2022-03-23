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


/// Algoritmos para obtener distancias y func básicas auxiliares

void dist_mat_op(float** coord, int num_points, float matrix[num_points][num_points]){
    for (int i=0; i<num_points;i++){
        matrix[i][i] = 0;
        for (int j=i; j<num_points;j++){
            matrix[i][j] = eucDisPP(coord[i],coord[j]);
            matrix[j][i] = matrix[i][j];
        }
    }
    return;
}

void swap_f(float *xp, float *yp){
    float temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap_i(int *xp, int *yp){
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}



/// Obtención de KNN y contrucción del grafo
void find_kmin_index(int n, float arr[n], int k, int k_min[k]){
    // Adaptación de l algoritmo de selection sort, encuentra los k valores minimos de un vector
    // y retorna sus índices.
    int i, j, min_idx;
    int idx_vec[n];
    float temp[n];

    for (i=0;i<n;i++){
        temp[i]=arr[i];
        idx_vec[i]=i;
    }

    for (i = 0; i < k; i++)
    {
        min_idx = i;
        for (j = i+1; j < n; j++)
          if (temp[j] < temp[min_idx])
            min_idx = j;
        swap_f(&temp[min_idx],&temp[i]);
        swap_i(&idx_vec[min_idx],&idx_vec[i]);
        k_min[i]=idx_vec[i];
    }
    return;
}

void knn_list(int n, float d_mat[n][n],  int k, int knn_mat[n][k]){
    for (int i=0; i<n;i++){
        find_kmin_index(n, d_mat[i], k, knn_mat[i]);
    }
    return;
}

bool n_in_array(int n, int dim ,int arr[dim]){
    int i=0;
    while(i<dim){
        if (n==arr[i]){
            return true;
        }
        i++;
    }
    return false;
}

void graph_mat(int n, float d_mat[n][n], int k, float graph[n][n]){
    int k_mat[n][k];
    knn_list(n,d_mat,k,k_mat);

    for (int row=0; row<n ; row++){
        for(int col=0; col<n ; col++)
            graph[row][col] = inf;
    }
    for (int row=0; row<n; row++){
        for(int col=0; col<n ; col++){
            if (n_in_array(col,k,k_mat[row])){
                graph[row][col] = d_mat[row][col];
                graph[col][row] = d_mat[row][col];
            }
        }
    }
    return;
}


/// Obtención de ruta más corta (Algoritmo de dijkstra)

int find_min_index(int n, float arr[n]){
    float min = arr[0];
    int min_idx = 0;
    for (int i = 1; i < n; i++)
    {
        if (min > arr[i]){
            min = arr[i];
            min_idx = i;
        }
    }
    return min_idx;
}

bool all_visited(int n, bool visited[n]){
    for ( int i=0; i<n;i++){
        if (!visited[i])
            return false;
    }
    return true;
}

bool all_infinites(int n, float cand[n]){
    for(int i=0;i<n;i++){
        if (cand[i] != inf)
            return false;
    }
    return true;
}

bool dijkstra( int n, float graph[n][n], int start, float short_path[n]){
    for (int i=0; i<n; i++){
        short_path[i] = inf;
    }
    bool visited[n];
    for (int i=0; i<n; i++){
        visited[i] = false;
    }
    short_path[start] = 0;

    while (!all_visited(n,visited)){
        float candidates[n];
        for (int i=0; i<n; i++){
            if (!visited[i])
                candidates[i] = short_path[i];
            else
                candidates[i] = inf;
        }

        if (all_infinites(n,candidates))
            return true;
        int m_idx = find_min_index(n,candidates);

        for (int i=0; i<n; i++){
            float new_dist = candidates[m_idx] + graph[m_idx][i];
            if (new_dist<short_path[i])
                short_path[i]=new_dist;
        }
        visited[m_idx] = true;

    }
    return false;
}

bool geodesic_mat(int n, float graph[n][n], float geo_dist[n][n]){
    for (int i = 0; i<n; i++){
         if(dijkstra(n,graph,i,geo_dist[i]))
            return true;
    }
    return false;
}

void fiber2geomat(int n, float** coord, int k, float geomat[n][n]){
    float dist[n][n];
    float graph[n][n];
    int knn = k+1;
    bool error = false;

    dist_mat_op(coord,n,dist);
    do{
        graph_mat(n,dist,knn,graph);
        error = geodesic_mat(n,graph,geomat);

        knn++;
    }while(error);

    return;
}


/// Aplicación del algoritmo de MDS

void centering_matrix(int n, float c_mat[n][n]){
    float trian = -1.0/n;
    float diag = (n-1.0)/n;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if (i==j)
                c_mat[i][j]=diag;
            else
                c_mat[i][j]=trian;
        }
    }

    return;
}

void sqr_dot_product(int dim, float M[dim][dim], float N[dim][dim],float X[dim][dim]){
    for (int i=0; i<dim; i++){
        for( int j=0; j<dim;j++){
            X[i][j]=0;
            for(int k=0;k<dim;k++){
                X[i][j] += M[i][k]*N[k][j];
            }
        }
    }
    return;
}

void sqr_matrix(int dim, float M[dim][dim]){
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            M[i][j] *= M[i][j];
        }
    }
}

void sqr_matrix_op(int dim, float M[dim][dim]){
    for(int i=0;i<dim;i++){
        for(int j=i+1;j<dim;j++){
            M[i][j] *= M[i][j];
            M[j][i] = M[i][j];
        }
    }
}

void dbl_centered_matrix( int n, float mat[n][n], float B[n][n]){
    sqr_matrix(n,mat);

    float C[n][n];
    centering_matrix(n,C);

    float C2[n][n];
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            C2[i][j] = (-0.5)*C[i][j];
    }

    float mul1[n][n];
    sqr_dot_product(n,C2,mat,mul1);
    sqr_dot_product(n,mul1,C,B);

    return;
}

void dbl_centered_matrix_op( int n, float mat[n][n], float B[n][n]){
    sqr_matrix_op(n,mat);

    float mean_col[n];
    float mean_row[n];

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            B[i][j]=mat[i][j];
        }
    }

    for(int i=0;i<n;i++){
        mean_col[i]=0;
        for(int j=0;j<n;j++){
            mean_col[i]+=B[j][i];
        }
        mean_col[i]/=n;
        for(int j=0;j<n;j++){
            B[j][i]-=mean_col[i];
        }
    }
    for(int i=0;i<n;i++){
        mean_row[i]=0;
        for(int j=0;j<n;j++){
            mean_row[i]+=B[i][j];
        }
        mean_row[i]/=n;
        for(int j=0;j<n;j++){
            B[i][j]-=mean_row[i];
        }
    }

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            B[i][j] *= (-0.5);
    }
    return;
}

void MDS_2D(int n, float** coord, int iter, float MDS[2][n], void(*dc_mat)(int,float[n][n],float[n][n])){
    float distmat[n][n];
    dist_mat_op(coord,n,distmat);

    float dcMat[n][n];
    dc_mat(n,distmat,dcMat);
    float eig_vec[n][n];
    float eig_val[n];

    int iter_num;
    int rot_num;

    jacobi_eigenvalue(n,dcMat,iter, eig_vec, eig_val,&iter_num, &rot_num);
    for(int i=0;i<n;i++){
        MDS[0][i] = eig_vec[0][i]*sqrt(eig_val[0]);
        MDS[1][i] = eig_vec[1][i]*sqrt(eig_val[1]);
    }
    return;
}

void ISOMAP_2D(int n, float** coord, int k, int iter, float isomap[2][n], void(*dc_mat)(int,float[n][n],float[n][n])){
    float geomat[n][n];
    fiber2geomat(n,coord,k,geomat);

    float dcMat[n][n];
    dc_mat(n,geomat,dcMat);

    float eig_vec[n][n];
    float eig_val[n];
    int iter_num;
    int rot_num;

    jacobi_eigenvalue(n, dcMat, iter, eig_vec, eig_val, &iter_num, &rot_num);

    for(int i=0;i<n;i++){
        isomap[0][i] = eig_vec[0][i]*sqrt(eig_val[0]);
        isomap[1][i] = eig_vec[1][i]*sqrt(eig_val[1]);
    }
    return;
}

float rad2deg(float rad){
    return rad*180.0/M_PI;
}

float deg2rad(float deg){
    return deg*M_PI/180.0;
}

void rotate_points(float rad, int n, float arr_2d[2][n]){
    float x,y;
    for(int i=0;i<n;i++){
         x = arr_2d[0][i]*cos(rad) - arr_2d[1][i]*sin(rad);
         y = arr_2d[0][i]*sin(rad) + arr_2d[1][i]*cos(rad);
         arr_2d[0][i]=x;
         arr_2d[1][i]=y;
    }
}

float rad_m(float p1[2], float p2[2]){
    float m = (p2[1]-p1[1])/(p2[0]-p1[0]);

    return atan(m);
};

void flip_ud(int n, float arr_2d[2][n]){
    for(int i=0;i<n;i++){
         arr_2d[1][i]= -arr_2d[1][i];
    }
}

void flip_lr(int n, float arr_2d[2][n]){
    for(int i=0;i<n;i++){
         arr_2d[0][i]= -arr_2d[1][i];
    }
}

/// valores propios (metodo de Jacobi)

/******************************************************************************/

void jacobi_eigenvalue ( int n, float a[n][n], int it_max, float v[n][n], float d[n], int *it_num, int *rot_num ){

/******************************************************************************/
/**
  Purpose:

    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.

  Discussion:

    This function computes the eigenvalues and eigenvectors of a
    real symmetric matrix, using Rutishauser's modfications of the classical
    Jacobi rotation method with threshold pivoting.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    C version by John Burkardt

    Input Parameters:

      int N, the order of the matrix.

      double A[N*N], the matrix, which must be square, real, and symmetric.

      int IT_MAX, the maximum number of iterations.


    Output Parameters:

      double V[N*N]: the matrix of eigenvectors.

      double D[N]: the eigenvalues, in descending order.

      int *IT_NUM: the total number of iterations.

      int *ROT_NUM: the total number of rotations.
*/

    float *bw;
    float c;
    float g;
    float gapq;
    float h;
    int i;
    int j;
    int k;
    int l;
    int m;
    int p;
    int q;
    float s;
    float t;
    float tau;
    float term;
    float termp;
    float termq;
    float theta;
    float thresh;
    float w;
    float *zw;

    r8mat_identity ( n, v );
    r8mat_diag_get_vector ( n, a, d );

    bw = ( float* ) malloc ( n * sizeof ( float) );
    zw = ( float* ) malloc ( n * sizeof ( float) );

    for ( i = 0; i < n; i++ ){
        bw[i] = d[i];
        zw[i] = 0.0;
    }
    *it_num = 0;
    *rot_num = 0;

    while ( *it_num < it_max ){
        *it_num = *it_num + 1;
        /*
          The convergence threshold is based on the size of the elements in
          the strict upper triangle of the matrix.
        */
        thresh = 0.0;
        for ( j = 0; j < n; j++ ){
            for ( i = 0; i < j; i++ ){
                thresh = thresh + a[j][i] * a[j][i] ;
            }
        }

        thresh = sqrt ( thresh ) / ( float) ( 4 * n );

        if ( thresh == 0.0 ){
            break;
        }

        for ( p = 0; p < n; p++ ){
            for ( q = p + 1; q < n; q++ ){
                gapq = 10.0 * fabs ( a[q][p] );
                termp = gapq + fabs ( d[p] );
                termq = gapq + fabs ( d[q] );

                //Annihilate tiny offdiagonal elements.
                if ( 4 < *it_num  &&  termp == fabs ( d[p] )  &&  termq == fabs ( d[q] ) ){
                    a[q][p] = 0.0;
                }


                //Otherwise, apply a rotation.
                else if ( thresh <= fabs ( a[q][p] ) ){
                    h = d[q] - d[p];
                    term = fabs ( h ) + gapq;

                    if ( term == fabs ( h ) ){
                        t = a[q][p] / h;
                    }
                    else{
                        theta = 0.5 * h / a[q][p];
                        t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
                        if ( theta < 0.0 ){
                            t = - t;
                        }
                    }
                    c = 1.0 / sqrt ( 1.0 + t * t );
                    s = t * c;
                    tau = s / ( 1.0 + c );
                    h = t * a[q][p];

                    //Accumulate corrections to diagonal elements.
                    zw[p] = zw[p] - h;
                    zw[q] = zw[q] + h;
                    d[p] = d[p] - h;
                    d[q] = d[q] + h;

                    a[q][p] = 0.0;

                    //Rotate, using information from the upper triangle of A only.
                    for ( j = 0; j < p; j++ ){
                        g = a[p][j];
                        h = a[q][j];
                        a[p][j] = g - s * ( h + g * tau );
                        a[q][j] = h + s * ( g - h * tau );
                    }

                    for ( j = p + 1; j < q; j++ ){
                        g = a[j][p];
                        h = a[q][j];
                        a[j][p] = g - s * ( h + g * tau );
                        a[q][j] = h + s * ( g - h * tau );
                    }

                    for ( j = q + 1; j < n; j++ ){
                        g = a[j][p];
                        h = a[j][q];
                        a[j][p] = g - s * ( h + g * tau );
                        a[j][q] = h + s * ( g - h * tau );
                    }

                    //Accumulate information in the eigenvector matrix.
                    for ( j = 0; j < n; j++ ){
                        g = v[p][j];
                        h = v[q][j];
                        v[p][j] = g - s * ( h + g * tau );
                        v[q][j] = h + s * ( g - h * tau );
                    }
                    *rot_num = *rot_num + 1;
                }
            }
        }

        for ( i = 0; i < n; i++ ){
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }

    //Restore upper triangle of input matrix.
    for ( j = 0; j < n; j++ ){
        for ( i = 0; i < j; i++ ){
          a[j][i] = a[i][j];
        }
    }

    //Descending sort the eigenvalues and eigenvectors.
    for ( k = 0; k < n - 1; k++ ){
        m = k;
        for ( l = k + 1; l < n; l++ ){
            if ( d[l] > d[m] ){
                m = l;
            }
        }

        if ( m != k ){
            t = d[m];
            d[m] = d[k];
            d[k] = t;
            for ( i = 0; i < n; i++ ){
                w = v[m][i];
                v[m][i] = v[k][i];
                v[k][i] = w;
            }
        }
    }

    free ( bw );
    free ( zw );

    return;
    /******************************************************************************/
}


void r8mat_diag_get_vector ( int n, float a[n][n],float v[n]){

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input, double A[N*N], the N by N matrix.

    Output, double V[N], the diagonal entries
    of the matrix.
*/

  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i][i];
  }

  return;
  /******************************************************************************/
}


void r8mat_identity  ( int n, float a[n][n] ){

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY sets an R8MAT to the identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/

    int i;
    int j;

    for ( j = 0; j < n; j++ ){
        for ( i = 0; i < n; i++ ){
            if ( i == j ){
                a[j][i] = 1.0;
            }
            else{
                a[j][i] = 0.0;
            }
        }
    }


    return;
    /******************************************************************************/
}


double r8mat_is_eigen_right ( int n, int k, double a[], double x[], double lambda[] ){

/******************************************************************************/
/*
  Purpose:

    R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.

  Discussion:

    An R8MAT is a matrix of doubles.

    This routine computes the Frobenius norm of

      A * X - X * LAMBDA

    where

      A is an N by N matrix,
      X is an N by K matrix (each of K columns is an eigenvector)
      LAMBDA is a K by K diagonal matrix of eigenvalues.

    This routine assumes that A, X and LAMBDA are all real.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, int K, the number of eigenvectors.
    K is usually 1 or N.

    Input, double A[N*N], the matrix.

    Input, double X[N*K], the K eigenvectors.

    Input, double LAMBDA[K], the K eigenvalues.

    Output, double R8MAT_IS_EIGEN_RIGHT, the Frobenius norm
    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
*/

  double *c;
  double error_frobenius;
  int i;
  int j;
  int l;

  c = ( double * ) malloc ( n * k * sizeof ( double ) );

  for ( j = 0; j < k; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( l = 0; l < n; l++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+l*n] * x[l+j*n];
      }
    }
  }

  for ( j = 0; j < k; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] - lambda[j] * x[i+j*n];
    }
  }

  error_frobenius = r8mat_norm_fro ( n, k, c );

  free ( c );

  return error_frobenius;
  /******************************************************************************/
}


double r8mat_norm_fro ( int m, int n, double a[] ){

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose Frobenius
    norm is desired.

    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
*/

  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
  /******************************************************************************/
}


void r8mat_print ( int m, int n, double a[], char *title ){

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/

  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
  /******************************************************************************/
}


void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title ){

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/

# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
/******************************************************************************/
}


void r8vec_print ( int n, double a[], char *title ){

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/

  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
  /******************************************************************************/
}


void timestamp ( ){

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/

# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


