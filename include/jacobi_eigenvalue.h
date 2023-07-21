void jacobi_eigenvalue ( int n, float a[n][n], int it_max, float v[n][n],
  float d[n], int *it_num, int *rot_num );

void r8mat_diag_get_vector ( int n, float a[n][n], float v[n] );

void r8mat_identity ( int n, float a[n][n] );

float r8mat_is_eigen_right ( int n, int k, float a[], float x[],
  float lambda[] );

float r8mat_norm_fro ( int m, int n, float a[] );

void r8mat_print ( int m, int n, float a[], char *title );

void r8mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, char *title );

void r8vec_print ( int n, float a[], char *title );

void timestamp ( void );

