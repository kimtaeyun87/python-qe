#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

// void invert_mat33(const double a[3][3], double b[3][3]);


double modulo(double x, double y)
{
    double z = fmod(x, y);
    // 1. z should be [0, y)
    if (z < 0) {z += y;}
    // 2. z ~ y --> z = 0
    if (fabs(z - y) < 1e-5) {z = 0.0;}
    return z;
}


double* invert_mat33(const double* A)
{
    double det = 0.0;

    det += A[0]*A[4]*A[8];
    det += A[3]*A[7]*A[2];
    det += A[6]*A[1]*A[5];
    det -= A[2]*A[4]*A[6];
    det -= A[1]*A[3]*A[8];
    det -= A[0]*A[7]*A[5];

    double* inv = (double*) malloc(sizeof(double)*9);

    inv[0] = (A[4]*A[8] - A[5]*A[7])/det;
    inv[1] = (A[2]*A[7] - A[8]*A[1])/det;
    inv[2] = (A[1]*A[5] - A[4]*A[2])/det;

    inv[3] = (A[5]*A[6] - A[8]*A[3])/det;
    inv[4] = (A[0]*A[8] - A[6]*A[2])/det;
    inv[5] = (A[2]*A[3] - A[5]*A[0])/det;

    inv[6] = (A[3]*A[7] - A[6]*A[4])/det;
    inv[7] = (A[1]*A[6] - A[7]*A[0])/det;
    inv[8] = (A[0]*A[4] - A[3]*A[1])/det;

    return inv;
}


int** get_equiv_atom_map(int M, const double* w, int N, const int* R, const double* t, int i0)
{
    int* map_I0 = (int*) malloc(sizeof(int)*M*N);
    int* map_I0_count = (int*) malloc(sizeof(int)*M);

    int max_count = 0;
    for (int m = 0; m < M; m++)
    {
        int count = 0;
        for (int n = 0; n < N; n++)
        {
            double w_m[3];
            for (int i = 0; i < 3; i++)
            {
                w_m[i] = t[3*n+i];
                for (int j = 0; j < 3; j++)
                {
                    w_m[i] += ((double) R[9*n+3*i+j])*w[3*i0+j];
                }
                w_m[i] = modulo(w_m[i], 1.0);
            }

            double y[3];
            for (int i = 0; i < 3; i++)
            {
                y[i] = w[3*m+i] - w_m[i];
            }

            double ynorm = cblas_dnrm2(3, y, 1);
            if (ynorm < 1e-5)
            {
                map_I0[N*m+count] = n;
                count += 1;
            }
        }

        if (count > max_count)
        {
            max_count = count;
        }
        map_I0_count[m] = count;
    }

    int* map_IJ = (int*) malloc(sizeof(int)*M*M*max_count);
    // int* map_IJ_count = (int*) malloc(sizeof(int)*M*M);

    for (int m1 = 0; m1 < M; m1++)
    {
        for (int m2 = 0; m2 < M; m2++)
        {
            // map_IJ_count[m1*M+m2] = 0;
            for (int count = 0; count < map_I0_count[m1]; count++)
            {
                int n = map_I0[N*m1+count];

                double* A = (double*) malloc(sizeof(double)*9);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        A[3*i+j] = (double) R[9*n+3*i+j];
                    }
                }
                double* A_inv = invert_mat33(A);

                double w_m3[3];
                for (int i = 0; i < 3; i++)
                {
                    w_m3[i] = 0.0;
                    for (int j = 0; j < 3; j++)
                    {
                        w_m3[i] += A_inv[3*i+j]*(w[3*m2+j] - t[3*n+j]);
                    }
                    w_m3[i] = modulo(w_m3[i], 1.0);
                }

                free(A);
                free(A_inv);

                for (int m3 = 0; m3 < M; m3++)
                {
                    double y[3];
                    for (int i = 0; i < 3; i++)
                    {
                        y[i] = w_m3[i] - w[3*m3+i];
                    }

                    double ynorm = cblas_dnrm2(3, y, 1);
                    if (ynorm < 1e-5)
                    {
                        map_IJ[max_count*M*m1+max_count*m2+count] = m3;
                        // map_IJ_count[M*m1+m2] += 1;
                    }
                }
            }
        }
    }
    
    int** map = (int**) malloc(sizeof(int*)*2);

    map[0] = map_I0_count;
    map[1] = map_IJ;

    free(map_I0);
    return map;
}

