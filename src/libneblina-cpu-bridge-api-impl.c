#ifndef __NEBLINAVECTOROCL
#define __NEBLINAVECTOROCL
#define BLOCK_DIM 16

#include <math.h>
#include "libneblina.h"
#include "neblina.h"
#include "bridge_api.h"

void InitEngine(int device)
{
    // printf("InitEngine\n");
    // InitCLEngine(device);
    // printf("end InitEngine\n");
}

void StopEngine()
{
    // ReleaseCLInfo(clinfo);
}

long get_Engine_Max_Memory_Allocation()
{
    return 0L;
}

void luDecomp(void *vector_a_dev, int n)
{

    return (void *)NULL;
}

double *addVectorF(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        // printf("%f %f\n",vector_a[i] , vector_b[i]);
        out[i] = vector_a[i] + vector_b[i];
    }

    return out;
}

void *addVectorC(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = vector_a[2 * i] + vector_b[2 * i];
        out[2 * i + 1] = vector_a[2 * i + 1] + vector_b[2 * i + 1];
    }

    return out;
}

void *addVectorFC(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = vector_a[i] + vector_b[2 * i];
        out[2 * i + 1] = vector_b[2 * i + 1];
    }

    return out;
}

void *prodVector(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[i] = vector_a[i] * vector_b[i];
    }

    return out;
}

void *vecConj(void *vector_a_dev, int n)
{
    return (void *)NULL;
}

void *vecConjugate(double *vector_a, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = vector_a[2 * i];
        out[2 * i + 1] = -vector_a[2 * i + 1];
    }
    return out;
}

void *vecAddOff2(void *vector_a_dev, int n)
{
    return (void *)NULL;
}

void *vecAddOff(double *vector_a, int offset, int parts)
{
    size_t n = parts * offset;

    double *out = (double *)malloc(offset * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < offset; i++)
    {
        double s = 0;
        for (int l = 0; l < parts; l++)
        {
            int idx = i + l * offset;
            s += vector_a[idx];
        }
        out[i] = s;
    }
    return out;
}

void *prodComplexVector(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        int idx_re = 2 * i;
        int idx_im = 2 * i + 1;
        out[idx_re] = vector_a[idx_re] * vector_b[idx_re] - vector_a[idx_im] * vector_b[idx_im];
        out[idx_im] = vector_a[idx_re] * vector_b[idx_im] + vector_a[idx_im] * vector_b[idx_re];
    }

    return out;
}

void *subVector(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        // printf("%f %f\n",vector_a[i] , vector_b[i]);
        out[i] = vector_a[i] - vector_b[i];
    }

    return out;
}

void *subVectorC(double *vector_a, double *vector_b, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = vector_a[2 * i] - vector_b[2 * i];
        out[2 * i + 1] = vector_a[2 * i + 1] - vector_b[2 * i + 1];
    }

    return out;
}

void *mulScalarVector(double *vector_a, double scalar, int n)
{
    double *out = (double *)malloc(n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[i] = scalar * vector_a[i];
    }

    return out;
}

void *mulComplexScalarVector(double *vector_a, double real, double imaginary, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = real * vector_a[i];
        out[2 * i + 1] = imaginary;
    }

    return out;
}

void *mulComplexScalarComplexVector(double *vector_a, double real, double imaginary, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = real * vector_a[2 * i];
        out[2 * i + 1] = imaginary * vector_a[2 * i + 1];
    }

    return out;
}

void *mulFloatScalarComplexVector(double *vector_a, double real, int n)
{
    double *out = (double *)malloc(2 * n * sizeof(double));
    // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        out[2 * i] = real * vector_a[2 * i];
        out[2 * i + 1] = vector_a[2 * i + 1];
    }

    return out;
}

void mulScalarMatRow(void *m, double scalar, int nrow, int ncols, int row)
{
    return (void *)NULL;
}

void mulScalarMatCol(void *m, double scalar, int nrow, int ncols, int col)
{
    return (void *)NULL;
}

void *matVecMul1(void *mDev, void *vDev, int ncols, int nrows)
{
    return (void *)NULL;
}

void *matVecMul2(void *mDev, void *vDev, int ncols, int nrows)
{
    return (void *)NULL;
}

int matrix_get_complex_real_index(int ncol, int i, int j)
{
    return 2 * (i * ncol + j);
}

int matrix_get_real_index(int ncol, int i, int j)
{
    return (i * ncol + j);
}

int matrix_get_complex_imag_index(int ncol, int i, int j)
{
    return matrix_get_complex_real_index(ncol, i, j) + 1;
}

void *matMulFloat(double *matrix_a, double *matrix_b, int nrows, int ncols, int ncol_matrix_a)
{
    double *out = (double *)malloc(nrows * ncols * sizeof(double));
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            double sum = 0;
            double vector_a;
            double vector_b;
#pragma omp unroll
            for (int k = 0; k < ncol_matrix_a; k++)
            {
                vector_a = matrix_a[i * ncol_matrix_a + k]; // matrix_a row
                vector_b = matrix_b[k * ncols + j];   // matrix_b col
                sum += vector_a * vector_b;
            }
            out[i * ncols + j] = sum;
        }
    }

    return out;
}

void *matMulComplex(double *matrix_a, double *matrix_b, int nrows, int ncols, int ncol_matrix_a)
{
    double *out = (double *)malloc(2 * nrows * ncols * sizeof(double));
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            int k;
            double sumre = 0, sumim = 0;
            double re1, im1, re2, im2;
#pragma omp unroll
            for (k = 0; k < ncol_matrix_a; k++)
            {
                int idx = matrix_get_complex_real_index(ncols, i, k);
                re1 = matrix_a[idx];
                im1 = matrix_a[idx + 1];

                idx = matrix_get_complex_real_index(ncols, k, j);
                re2 = matrix_b[idx];
                im2 = matrix_b[idx + 1];

                sumre += re1 * re2 - im1 * im2;
                sumim += re1 * im2 + re2 * im1;
            }
            int idx_out = matrix_get_complex_real_index(ncols, i, j);
            out[idx_out] = sumre;
            out[idx_out + 1] = sumim;
        }
    }

    return out;
}

void *matMulFloatComplex(double *matrix_a, double *matrix_b, int nrows, int ncols, int ncol_matrix_a)
{
    double *out = (double *)malloc(2 * nrows * ncols * sizeof(double));
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            int k;
            double sumre = 0, sumim = 0;
            double re1, re2, im2;
#pragma omp unroll
            for (k = 0; k < ncol_matrix_a; k++)
            {
                int idx = matrix_get_complex_real_index(ncol_matrix_a, i, k);
                re1 = matrix_a[idx];

                idx = matrix_get_complex_real_index(ncols, k, j);
                re2 = matrix_b[idx];
                im2 = matrix_b[idx + 1];

                sumre += re1 * re2;
                sumim += re1 * im2;
            }
            int idx_out = matrix_get_complex_real_index(ncols, i, j);
            out[idx_out] = sumre;
            out[idx_out + 1] = sumim;
        }
    }
    return out;
}

void *matMul(double *matrix_a, double *matrix_b, int nrows, int ncols, int ncol_matrix_a, int atype, int btype)
{
    if (atype == T_FLOAT && btype == T_FLOAT)
    {
        return matMulFloat(matrix_a, matrix_b, nrows, ncols, ncol_matrix_a);
    }
    else if (atype == T_COMPLEX && btype == T_COMPLEX)
    {
        return matMulComplex(matrix_a, matrix_b, nrows, ncols, ncol_matrix_a);
    }
    else if (atype == T_FLOAT && btype == T_COMPLEX)
    {
        return matMulFloatComplex(matrix_a, matrix_b, nrows, ncols, ncol_matrix_a);
    }
}

void *matMul2(void *m1Dev, void *m2Dev, int nrows, int ncols, int qq)
{
    return (void *)NULL;
}

void matSquare(void **outLin, void **idxOutLin,
               void **outCol, void **idxOutCol,
               void *mLin, void *idxLin,
               void *mCol, void *idxCol,
               int maxcols, int N)
{
    return (void *)NULL;
}

void *matVecMul3(double *matrix, double *vector, int ncols, int nrows)
{
    double *out = (double *)malloc(nrows * sizeof(double));

#pragma omp parallel for
    for (int i = 0; i < nrows; i++)
    {
        double sum = 0;
#pragma omp unroll
        for (int j = 0; j < ncols; j++)
        {
            double vector_a;
            double vector_b;
            int idx1 = matrix_get_real_index(ncols, i, j);
            vector_a = matrix[idx1];
            vector_b = vector[j];
            sum += vector_a * vector_b;
        }
        int idx_out = i; // matrix_get_real_index(ncols, i, j);
        out[idx_out] = sum;
    }
    return out;
}

void *sparseVecMul(void *mDev, void *idxCol, void *vDev, int nrows, int maxCols)
{
    double *vec_out = (double *)malloc(nrows * maxCols * sizeof(double));
    double sum_re = 0, re_m, re_v;
    double *m = (double *)mDev;
    double *vec_in = (double *)vDev;
    int *col_idx = (int *)idxCol;

#pragma omp parallel for
    for (int idx = 0; idx < nrows; idx++)
    {
        double sum = 0;
        int row = idx;
        int midx = idx * maxCols; // midx -> M de max - m_idx ficaria melhor
        int i;
        for (i = 0; i < maxCols; i++)
        {
            int col = col_idx[midx + i];
            sum += (col != -1) ? m[midx + i] * vec_in[col] : 0; // -1 para pular de linha
        }
        vec_out[row] = sum;
    }

    return (void *)vec_out;
}
void *sparseComplexVecMul(void *mDev, void *idxCol, void *vDev, int nrows, int maxCols)
{
    double *vec_out = (double *)malloc(2 * nrows * maxCols * sizeof(double));

#pragma omp parallel for
    for (size_t idx = 0; idx < nrows; idx++)
    {
        double sum_re = 0, sum_im = 0, re_m, im_m, re_v, im_v;
        double *m = (double *)mDev;
        double *vec_in = (double *)vDev;
        int *col_idx = (int *)idxCol;
        int row = idx, i, col, idxt;
        for (i = 0; i < maxCols; i++)
        {
            idxt = (row * maxCols) + i;
            col = col_idx[idxt];
            if (col == -1)
                continue;
            re_m = m[2 * idxt];
            im_m = m[2 * idxt + 1];
            re_v = vec_in[2 * col];
            im_v = vec_in[2 * col + 1];
            sum_re += re_m * re_v - im_m * im_v;
            sum_im += re_m * im_v + im_m * re_v;
        }
        vec_out[2 * row] = sum_re;
        vec_out[2 * row + 1] = sum_im;
    }
    return (void *)vec_out;
}

void *matVecMul3Complex(double *matrix, double *vector, int ncols, int nrows)
{
    double *out = (double *)malloc(2 * nrows * sizeof(double));

#pragma omp parallel for
    for (int i = 0; i < nrows; i++)
    {
        double sumre = 0, sumim = 0;
        for (int j = 0; j < ncols; j++)
        {

            double re1, im1, re2, im2;
            int idx = matrix_get_complex_real_index(ncols, i, j);
            re1 = matrix[idx];
            im1 = matrix[idx + 1];

            idx = 2 * j;
            re2 = vector[idx];
            im2 = vector[idx + 1];

            sumre += (re1 * re2) - (im1 * im2);
            sumim += (re1 * im2) + (re2 * im1);
        }
        int idx_out = 2 * i; // matrix_get_complex_real_index(ncols,i,j);
        out[idx_out] = sumre;
        out[idx_out + 1] = sumim;
    }

    return out;
}

void *matTranspose(void *mDev, int ncols, int nrows)
{
    return (void *)NULL;
}

double sumVector(double *vector_a, int len)
{
    double out = 0;
    for (int i = 0; i < len; i++)
    {
        out += vector_a[i];
    }

    return out;
}

double normVector(void *vDev, int len)
{
    return 0.0;
}

double dotVector(void *vector_a_dev, void *v2Dev, int len)
{
    double sum = 0;

    double *vector_a = (double *)vector_a_dev;
    double *vector_b = (double *)v2Dev;
    for (int i = 0; i < len; i++)
    {
        sum += vector_a[i] * vector_b[i];
    }

    return sum;
}

void dotVectorComplex(double *out_re, double *out_im, void *vector_a_dev, void *v2Dev, int len)
{
    return (void *)NULL;
}

#endif
