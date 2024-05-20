#include "libneblina.h"
#include <stdio.h>
#include <stdlib.h>

matrix_t *matrix_new(int nrow, int ncol, data_type type, int initialize, void *data)
{
    matrix_t *ret = (matrix_t *)malloc(sizeof(matrix_t));

    if (initialize && data == NULL)
    {
        if (type == T_INT)
        {
            ret->value.i = (int *)malloc(nrow * ncol * sizeof(int));
        }
        else if (type == T_FLOAT)
        {
            ret->value.f = (double *)malloc(nrow * ncol * sizeof(double));
        }
        else if (type == T_COMPLEX)
        {
            ret->value.f = (double *)malloc(2 * nrow * ncol * sizeof(double));
        }
        ret->externalData = 0;
    }
    else if (data != NULL)
    {
        ret->value.f = (double *)data;
        ret->externalData = 1;
    }
    else
    {
        ret->value.f = NULL;
        ret->externalData = 0;
    }

    ret->type = type;
    ret->nrow = nrow;
    ret->ncol = ncol;
    ret->location = LOCHOS;
    ret->extra = NULL;

    return ret;
}

void matrix_delete(matrix_t *matrix)
{
    if (matrix->value.f != NULL && matrix->externalData == 0)
    {
        free(matrix->value.f);
    }
    else if (matrix->extra != NULL && matrix->externalData == 0)
    {
        free(matrix->extra);
    }
    free(matrix);
}

void matreqhost(matrix_t *v)
{
    if (v->location == LOCHOS)
        return;

    v->location = LOCHOS;
    v->value.f = v->extra;
    v->extra = NULL;
}

void matreqdev(matrix_t *v)
{
    if (v->location == LOCDEV)
        return;

    v->location = LOCDEV;
    v->extra = v->value.f;
    v->value.f = NULL;
}
void matrix_set_real_value(matrix_t *matrix, int i, int j, double r)
{
    matrix->value.f[i * matrix->ncol + j] = r;
}

double matrix_get_real_value(matrix_t *matrix, int i, int j)
{
    return matrix->value.f[i * matrix->ncol + j];
}

void matrix_set_complex_value(matrix_t *matrix, int i, int j, double r, double im)
{
    int idx = 2 * (i * matrix->ncol + j);
    matrix->value.f[idx] = r;
    matrix->value.f[idx + 1] = im;
}

double matrix_get_complex_real_value(matrix_t *matrix, int i, int j)
{
    int idx = 2 * (i * matrix->ncol + j);
    return matrix->value.f[idx];
}

double matrix_get_complex_imaginary_value(matrix_t *matrix, int i, int j)
{
    int idx = 2 * (i * matrix->ncol + j);
    return matrix->value.f[idx + 1];
}

// if we want a part of the column we can simply state a size
// smaller than the size of the column
double *matrix_copy_col(matrix_t *matrix, int j, int ini, int size)
{
    if (size > matrix->ncol)
    {
        size = matrix->ncol;
    }
    double *col = (double *)malloc(size * sizeof(double));
    int idx = 0;
    for (int k = ini; k < (size + ini); k++)
    {
        col[idx] = matrix->value.f[k * matrix->ncol + j];
        idx++;
    }

    return col;
}

// if we want a part of the row we can simply state a size
// smaller than the size of the row
double *matrix_copy_row(matrix_t *matrix, int i, int ini, int size)
{
    if (size > matrix->nrow)
    {
        size = matrix->nrow;
    }
    double *row = (double *)malloc(size * sizeof(double));
    int idx = 0;

    for (int k = ini; k < (size + ini); k++)
    {
        row[idx] = matrix->value.f[i * matrix->ncol + k];
        idx++;
    }
    return row;
}