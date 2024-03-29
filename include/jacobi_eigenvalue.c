#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "jacobi_eigenvalue.h"

/******************************************************************************/

void jacobi_eigenvalue(int n, float a[n][n], int it_max, float v[n][n],
                       float d[n], int *it_num, int *rot_num)

/******************************************************************************/
/*
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

    John Burkardt

  Reference:

    Gene Golub, Charles VanLoan,
    Matrix Computations,
    Third Edition,
    Johns Hopkins, 1996,
    ISBN: 0-8018-4513-X,
    LC: QA188.G65.

  Input:

    int N, the order of the matrix.

    float A[N*N], the matrix, which must be square, real,
    and symmetric.

    int IT_MAX, the maximum number of iterations.

  Output:

    float V[N*N], the matrix of eigenvectors.

    float D[N], the eigenvalues, in descending order.

    int *IT_NUM, the total number of iterations.

    int *ROT_NUM, the total number of rotations.
*/
{
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

  r8mat_identity(n, v);

  r8mat_diag_get_vector(n, a, d);

  bw = (float *)malloc(n * sizeof(float));
  zw = (float *)malloc(n * sizeof(float));

  for (i = 0; i < n; i++)
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  *it_num = 0;
  *rot_num = 0;

  while (*it_num < it_max)
  {
    *it_num = *it_num + 1;
    /*
      The convergence threshold is based on the size of the elements in
      the strict upper triangle of the matrix.
    */
    thresh = 0.0;
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < j; i++)
      {
        thresh = thresh + a[j][i] * a[j][i];
      }
    }

    thresh = sqrt(thresh) / (float)(4 * n);

    if (thresh == 0.0)
    {
      break;
    }

    for (p = 0; p < n; p++)
    {
      for (q = p + 1; q < n; q++)
      {
        gapq = 10.0 * fabs(a[q][p]);
        termp = gapq + fabs(d[p]);
        termq = gapq + fabs(d[q]);
        /*
          Annihilate tiny offdiagonal elements.
        */
        if (4 < *it_num &&
            termp == fabs(d[p]) &&
            termq == fabs(d[q]))
        {
          a[q][p] = 0.0;
        }
        /*
          Otherwise, apply a rotation.
        */
        else if (thresh <= fabs(a[q][p]))
        {
          h = d[q] - d[p];
          term = fabs(h) + gapq;

          if (term == fabs(h))
          {
            t = a[q][p] / h;
          }
          else
          {
            theta = 0.5 * h / a[q][p];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
            {
              t = -t;
            }
          }
          c = 1.0 / sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[q][p];
          /*
            Accumulate corrections to diagonal elements.
          */
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[q][p] = 0.0;
          /*
            Rotate, using information from the upper triangle of A only.
          */
          for (j = 0; j < p; j++)
          {
            g = a[p][j];
            h = a[q][j];
            a[p][j] = g - s * (h + g * tau);
            a[q][j] = h + s * (g - h * tau);
          }

          for (j = p + 1; j < q; j++)
          {
            g = a[j][p];
            h = a[q][j];
            a[j][p] = g - s * (h + g * tau);
            a[q][j] = h + s * (g - h * tau);
          }

          for (j = q + 1; j < n; j++)
          {
            g = a[j][p];
            h = a[j][p];
            a[j][p] = g - s * (h + g * tau);
            a[j][p] = h + s * (g - h * tau);
          }
          /*
            Accumulate information in the eigenvector matrix.
          */
          for (j = 0; j < n; j++)
          {
            g = v[p][j];
            h = v[q][j];
            v[p][j] = g - s * (h + g * tau);
            v[q][j] = h + s * (g - h * tau);
          }
          *rot_num = *rot_num + 1;
        }
      }
    }

    for (i = 0; i < n; i++)
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
  /*
    Restore upper triangle of input matrix.
  */
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      a[j][i] = a[i][j];
    }
  }
  /*
    Descending sort the eigenvalues and eigenvectors.
  */
  for (k = 0; k < n - 1; k++)
  {
    m = k;
    for (l = k + 1; l < n; l++)
    {
      if (d[l] > d[m])
      {
        m = l;
      }
    }

    if (m != k)
    {
      t = d[m];
      d[m] = d[k];
      d[k] = t;
      for (i = 0; i < n; i++)
      {
        w = v[m][i];
        v[m][i] = v[k][i];
        v[k][i] = w;
      }
    }
  }

  free(bw);
  free(zw);

  return;
}
/******************************************************************************/

void r8mat_diag_get_vector(int n, float a[n][n], float v[n])

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

  Input:

    int N, the number of rows and columns of the matrix.

    float A[N*N], the N by N matrix.

  Output:

    float V[N], the diagonal entries
    of the matrix.
*/
{
  int i;

  for (i = 0; i < n; i++)
  {
    v[i] = a[i][i];
  }

  return;
}
/******************************************************************************/

void r8mat_identity(int n, float a[n][n])

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

  Input:

    int N, the order of A.

  Output:

    float A[N*N], the N by N identity matrix.
*/
{
  int i;
  int j;

  for (j = 0; j < n; j++)
  {
    a[j][j] = 1.0;
    for (i = j + 1; i < n; i++)
    {
      a[j][i] = 0.0;
      a[i][j] = 0.0;
    }
  }

  return;
}
/******************************************************************************/

float r8mat_is_eigen_right(int n, int k, float a[], float x[],
                           float lambda[])

/******************************************************************************/
/*
  Purpose:

    R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.

  Discussion:

    An R8MAT is a matrix of floats.

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

  Input:

    int N, the order of the matrix.

    int K, the number of eigenvectors.
    K is usually 1 or N.

    float A[N*N], the matrix.

    float X[N*K], the K eigenvectors.

    float LAMBDA[K], the K eigenvalues.

  Output:

    float R8MAT_IS_EIGEN_RIGHT, the Frobenius norm
    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
*/
{
  float *c;
  float error_frobenius;
  int i;
  int j;
  int l;

  c = (float *)malloc(n * k * sizeof(float));

  for (j = 0; j < k; j++)
  {
    for (i = 0; i < n; i++)
    {
      c[i + j * n] = 0.0;
      for (l = 0; l < n; l++)
      {
        c[i + j * n] = c[i + j * n] + a[i + l * n] * x[l + j * n];
      }
    }
  }

  for (j = 0; j < k; j++)
  {
    for (i = 0; i < n; i++)
    {
      c[i + j * n] = c[i + j * n] - lambda[j] * x[i + j * n];
    }
  }

  error_frobenius = r8mat_norm_fro(n, k, c);

  free(c);

  return error_frobenius;
}
/******************************************************************************/

float r8mat_norm_fro(int m, int n, float a[])

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

  Input:

    int M, the number of rows in A.

    int N, the number of columns in A.

    float A[M*N], the matrix whose Frobenius
    norm is desired.

  Output:

    float R8MAT_NORM_FRO, the Frobenius norm of A.
*/
{
  int i;
  int j;
  float value;

  value = 0.0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      value = value + pow(a[i + j * m], 2);
    }
  }
  value = sqrt(value);

  return value;
}
/******************************************************************************/

void r8mat_print(int m, int n, float a[], char *title)

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

  Input:

    int M, the number of rows in A.

    int N, the number of columns in A.

    float A[M*N], the M by N matrix.

    char *TITLE, a title.
*/
{
  r8mat_print_some(m, n, a, 1, 1, m, n, title);

  return;
}
/******************************************************************************/

void r8mat_print_some(int m, int n, float a[], int ilo, int jlo, int ihi,
                      int jhi, char *title)

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

  Input:

    int M, the number of rows of the matrix.
    M must be positive.

    int N, the number of columns of the matrix.
    N must be positive.

    float A[M*N], the matrix.

    int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    char *TITLE, a title.
*/
{
#define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf(stdout, "\n");
  fprintf(stdout, "%s\n", title);

  if (m <= 0 || n <= 0)
  {
    fprintf(stdout, "\n");
    fprintf(stdout, "  (None)\n");
    return;
  }
  /*
    Print the columns of the matrix, in strips of 5.
  */
  for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
  {
    j2hi = j2lo + INCX - 1;
    if (n < j2hi)
    {
      j2hi = n;
    }
    if (jhi < j2hi)
    {
      j2hi = jhi;
    }

    fprintf(stdout, "\n");
    /*
      For each column J in the current range...

      Write the header.
    */
    fprintf(stdout, "  Col:  ");
    for (j = j2lo; j <= j2hi; j++)
    {
      fprintf(stdout, "  %7d     ", j - 1);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "  Row\n");
    fprintf(stdout, "\n");
    /*
      Determine the range of the rows in this strip.
    */
    if (1 < ilo)
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if (m < ihi)
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for (i = i2lo; i <= i2hi; i++)
    {
      /*
        Print out (up to) 5 entries in row I, that lie in the current strip.
      */
      fprintf(stdout, "%5d:", i - 1);
      for (j = j2lo; j <= j2hi; j++)
      {
        fprintf(stdout, "  %14f", a[i - 1 + (j - 1) * m]);
      }
      fprintf(stdout, "\n");
    }
  }

  return;
#undef INCX
}
/******************************************************************************/

void r8vec_print(int n, float a[], char *title)

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

  Input:

    int N, the number of components of the vector.

    float A[N], the vector to be printed.

    char *TITLE, a title.
*/
{
  int i;

  fprintf(stdout, "\n");
  fprintf(stdout, "%s\n", title);
  fprintf(stdout, "\n");
  for (i = 0; i < n; i++)
  {
    fprintf(stdout, "  %8d: %14f\n", i, a[i]);
  }

  return;
}
/******************************************************************************/

void timestamp()

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
*/
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  fprintf(stdout, "%s\n", time_buffer);

  return;
#undef TIME_SIZE
}
