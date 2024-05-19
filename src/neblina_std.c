#include "neblina_std.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>

#include "libneblina.h"
#include "neblina.h"
#include "neblina_list.h"
#include "bridge_api.h"

#define TRUE 1

void runerror(char *strerr)
{
    fprintf(stderr, " runtime error: %s\n", strerr);
    exit(1);
}

void **mat_len_col(void **i, int *status)
{
    object_t out; // = (object_t *) malloc( sizeof( object_t ) );
    object_t **in = (object_t **)i;
    type(out) = T_INT;
    matrix_t *mat = (matrix_t *)vvalue(*in[0]);
    ivalue(out) = mat->ncol;
    static void *ret[1];
    clear_input(i, 1);
    ret[0] = (void *)&out;
    if (status != NULL)
    {
        *status = 0;
    }
    return ret;
}

void **mat_len_row(void **i, int *status)
{
    object_t **in = (object_t **)i;
    if (status != NULL)
    {
        *status = 0;
    }
    if (type(*in[0]) == T_MATRIX)
    {
        object_t out; // = (object_t *) malloc( sizeof( object_t ) );

        type(out) = T_INT;
        matrix_t *mat = (matrix_t *)vvalue(*in[0]);
        ivalue(out) = mat->nrow;
        static void *ret[1];
        clear_input(i, 1);
        ret[0] = (void *)&out;
        return ret;
    }
    else
    {
        if (status != NULL)
        {
            *status = -1;
        }
        runerror("Runtime error: no matrix found\n");
    }
    return NULL;
}

void delete_object_array(object_t **in, int len)
{
    if (in != NULL)
    {
        for (int i = 0; i < len; i++)
        {
            delete_object(in[i]);
        }
        free(in);
    }
}

void delete_object(object_t *in)
{
    if (in != NULL)
    {
        free(in);
    }
}

void **vec_add(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    vector_t *b = (vector_t *)vvalue(*in[1]);
    vector_t *r = m->bridges[index].vector_new(b->len, b->type, 0, NULL);
    // apenas para CPU
    //  free(r->value.f);
    m->bridges[index].vecreqdev(a);
    m->bridges[index].vecreqdev(b);
    m->bridges[index].vecreqdev(r);

    if (b->type == T_FLOAT)
    {
        r->extra = (void *)m->bridges[index].addVectorF_f(a->extra, b->extra, b->len);
    }
    else if (b->type == T_COMPLEX)
    {
        r->extra = (void *)m->bridges[index].addVectorC_f(a->extra, b->extra, b->len);
    }

    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **vec_conj(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    vector_t *r = m->bridges[index].vector_new(a->len, T_COMPLEX, 0, NULL);
    // apenas para cpu
    //  free( r->value.f);
    m->bridges[index].vecreqdev(a);
    m->bridges[index].vecreqdev(r);

    r->extra = (void *)m->bridges[index].vecConjugate_f(a->extra, a->len);

    clear_input(i, 1);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **vec_prod(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    vector_t *b = (vector_t *)vvalue(*in[1]);
    m->bridges[index].vecreqdev(a);
    m->bridges[index].vecreqdev(b);
    vector_t *r = (vector_t *)malloc(sizeof(vector_t));
    if (a->type == T_FLOAT)
    {
        r->extra = (void *)m->bridges[index].prodVector_f(a->extra, b->extra, b->len);
    }
    else
    {
        r->extra = (void *)m->bridges[index].prodComplexVector_f(a->extra, b->extra, b->len);
    }

    r->len = b->len;
    r->type = a->type;
    r->location = LOCDEV;
    r->value.f = NULL;
    clear_input(i, 2);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **vec_sum(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    object_t *out = (object_t *)malloc(sizeof(object_t));

    m->bridges[index].vecreqdev(a);
    out->value.f = m->bridges[index].sumVector_f(a->extra, a->len);
    out->type = T_FLOAT;

    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)out;
}

void **vec_norm(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    object_t *out = (object_t *)malloc(sizeof(object_t));
    m->bridges[index].vecreqdev(a);
    fvalue(*out) = m->bridges[index].normVector_f(a->extra, a->len);
    type(*out) = T_FLOAT;
    static void *ret[1];
    ret[0] = (void *)out;
    if (status != NULL)
    {
        *status = 0;
    }
    return ret;
}

void **vec_dot(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *v1 = (vector_t *)vvalue(*in[0]);
    vector_t *v2 = (vector_t *)vvalue(*in[1]);
    m->bridges[index].vecreqdev(v1);
    m->bridges[index].vecreqdev(v2);
    object_t *out;
    if (type(*v1) == T_FLOAT && type(*v2) == T_FLOAT)
    {
        out = (object_t *)malloc(sizeof(object_t));
        fvalue(*out) = m->bridges[index].dotVector_f(v1->extra, v2->extra, v1->len);
        type(*out) = T_FLOAT;
    }
    else
    {
        out = (object_t *)malloc(sizeof(object_t));
        double re, im;
        m->bridges[index].dotVectorComplex_f(&re, &im, v1->extra, v2->extra, v1->len);
        complex_t *res = (complex_t *)malloc(sizeof(complex_t));
        res->im = im;
        res->re = re;
        vvalue(*out) = (void *)res;
        type(*out) = T_COMPLEX;
    }

    static void *ret[1];
    ret[0] = (void *)out;
    if (status != NULL)
    {
        *status = 0;
    }
    return ret;
}

void **vec_sub(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *a = (vector_t *)vvalue(*in[0]);
    vector_t *b = (vector_t *)vvalue(*in[1]);
    vector_t *r = m->bridges[index].vector_new(b->len, b->type, 0, NULL);
    // apenas para cpu
    // free(r->value.f);
    m->bridges[index].vecreqdev(a);
    m->bridges[index].vecreqdev(b);
    m->bridges[index].vecreqdev(r);

    if (b->type == T_FLOAT)
    {
        r->extra = (void *)m->bridges[index].subVector_f(a->extra, b->extra, b->len);
    }
    else if (b->type == T_COMPLEX)
    {
        r->extra = (void *)m->bridges[index].subVectorC_f(a->extra, b->extra, b->len);
    }
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **mat_add(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    matrix_t *a = (matrix_t *)vvalue(*in[0]);
    matrix_t *b = (matrix_t *)vvalue(*in[1]);
    m->bridges[index].matreqdev(a);
    m->bridges[index].matreqdev(b);
    // object_t out;
    matrix_t *r;
    if (a->type == T_FLOAT && b->type == T_FLOAT)
    {
        r = m->bridges[index].matrix_new(b->ncol, b->nrow, T_FLOAT, 0, NULL);
        r->extra = m->bridges[index].addVectorF_f(a->extra, b->extra, b->nrow * b->ncol);
        r->location = LOCDEV;
    }
    else if (a->type == T_COMPLEX && b->type == T_COMPLEX)
    {
        r = m->bridges[index].matrix_new(b->ncol, b->nrow, T_COMPLEX, 0, NULL);
        r->extra = m->bridges[index].addVectorC_f(a->extra, b->extra, b->nrow * b->ncol);
        r->location = LOCDEV;
    }
    else if ((a->type == T_FLOAT && b->type == T_COMPLEX) ||
             (a->type == T_COMPLEX && b->type == T_FLOAT))
    {
        r = m->bridges[index].matrix_new(b->ncol, b->nrow, T_COMPLEX, 0, NULL);
        r->ncol = b->ncol;
        r->nrow = b->nrow;
        r->type = T_COMPLEX;
        if (a->type == T_FLOAT)
        {
            r->extra = m->bridges[index].addVectorFC_f(a->extra, b->extra, b->nrow * b->ncol);
        }
        else
        {
            r->extra = m->bridges[index].addVectorFC_f(b->extra, a->extra, b->nrow * b->ncol);
        }
        r->location = LOCDEV;
        r->value.f = NULL;
    }
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **mat_sub(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    matrix_t *a = (matrix_t *)vvalue(*in[1]);
    matrix_t *b = (matrix_t *)vvalue(*in[0]);
    object_t out; // = (object_t *) malloc( sizeof( object_t ) );
    matrix_t *r = (matrix_t *)malloc(sizeof(matrix_t));

    m->bridges[index].matreqdev(a);
    m->bridges[index].matreqdev(b);
    
    r->ncol = a->ncol;
    r->nrow = a->nrow;
    r->type = T_FLOAT;
    r->extra = m->bridges[index].subVector_f(a->extra, b->extra, b->nrow * b->ncol);
    r->location = LOCDEV;
    r->value.f = NULL;
    
    type(out) = T_MATRIX;
    vvalue(out) = (void *)r;
    static void *ret[1];
    ret[0] = (void *)&out;
    clear_input(i, 2);
    
    if (status != NULL)
    {
        *status = 0;
    }
    
    return ret;
}

void **mat_mul(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    matrix_t *a = (matrix_t *)vvalue(*in[0]);
    matrix_t *b = (matrix_t *)vvalue(*in[1]);

    matrix_t *r;
    long matrix_size = 0;
    if (a->type == T_FLOAT && b->type == T_FLOAT)
    {
        r = m->bridges[index].matrix_new(b->ncol, a->nrow, T_FLOAT, 0, NULL);
        r->type = T_FLOAT;
        matrix_size = b->ncol;
    }
    else if ((a->type == T_COMPLEX && b->type == T_COMPLEX) ||
             (a->type == T_FLOAT && b->type == T_COMPLEX))
    {
        r = m->bridges[index].matrix_new(b->ncol, a->nrow, T_COMPLEX, 0, NULL);
        r->type = T_COMPLEX;
        matrix_size = 2 * b->ncol;
    }
    else
    {
        runerror("Invalid types for mat_mul\n");
    }
    struct timeval stop, start, ini, end, tval_result;
    long max_mem = m->bridges[index].get_Engine_Max_Memory_Allocation_f();
    printf("matrix_size=%ld max_mem=%ld (max_mem / sizeof(double))=%ld\n", matrix_size, max_mem, (max_mem / sizeof(double)));
    printf("matrix_size < (max_mem / sizeof(double))=%d\n", matrix_size < (max_mem / sizeof(double)));
    if (TRUE)
    { //|| (matrix_size * matrix_size) < (max_mem / sizeof(double))
        // gettimeofday(&ini, NULL);
        m->bridges[index].matreqdev(a);
        m->bridges[index].matreqdev(b);

        r->extra = (void *)m->bridges[index].matMul_f(a->extra, b->extra, a->nrow, b->ncol, a->ncol, a->type, b->type);
        // gettimeofday(&end, NULL);
        // timersub(&end, &ini, &tval_result);
        // printf("Time elapsed: %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
        r->location = LOCDEV;
        r->value.f = NULL;
    }
    else
    {
        // A * B -> row from A and col from B
        // use max_mem
        printf("processing chunks\n");
        r->location = LOCHOS;
        r->value.f = (double *)calloc(a->nrow * b->ncol, sizeof(double));
        printf("callocated\n");

        printf("1\n");
        printf("max_mem=%d\n", max_mem);
        printf("(max_mem / sizeof(double))=%d\n", (max_mem / sizeof(double)));
        long qty_chunks = ceil((matrix_size * 1.0) / (max_mem / sizeof(double)));
        printf("qty_chunks=%ld\n", qty_chunks);
        printf("2\n");
        long chunk_size = matrix_size / qty_chunks;
        printf("3\n");
        printf("chunk_size=%ld\n", chunk_size);

        for (int j = 0; j < b->ncol; j++)
        { // we will move through columns first
            gettimeofday(&ini, NULL);
            for (int c = 0; c < qty_chunks; c++)
            { // then we move through the chunks
                // gettimeofday(&start, NULL);
                double *B_col = m->bridges[index].matrix_copy_col(b, j, c * chunk_size, chunk_size);
                // gettimeofday(&stop, NULL);
                // printf(" col copy took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
                // getchar();
                // printf("B_col=%p\n",B_col);
                // for (int n=0; n <chunk_size; n++){
                //     printf("B_col[%d]=%lf\n", n, B_col[n]);
                // }
                vector_t *col = m->bridges[index].vector_new(chunk_size, T_FLOAT, 0, B_col);
                // for (int n=0; n <col->len; n++){
                //     printf("col[%d]=%lf\n", n, col->value.f[n]);
                //     printf("B_col[%d]=%lf\n", n, B_col[n]);
                // }
                m->bridges[index].vecreqdev(col);
                for (int i = 0; i < a->nrow; i++)
                {
                    // for the same column chunk that was copied to device memory
                    // we calculate the dot product for all the row chunks to leverage
                    // the column that was copied first (columns will inccur in more
                    // cache misses on the CPU)

                    // gettimeofday(&start, NULL);
                    double *A_row = m->bridges[index].matrix_copy_row(a, i, c * chunk_size, chunk_size);
                    // gettimeofday(&stop, NULL);
                    // printf("  row copy took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
                    // printf("A_row=%p\n",A_row);
                    vector_t *row = m->bridges[index].vector_new(chunk_size, T_FLOAT, 0, A_row);
                    // for (int n=0; n <row->len; n++){
                    //     printf("row[%d]=%lf\n", n, row->value.f[n]);
                    //     printf("A_row[%d]=%lf\n", n, A_row[n]);
                    // }
                    // gettimeofday(&start, NULL);
                    m->bridges[index].vecreqdev(row);
                    // gettimeofday(&stop, NULL);
                    // printf("  row move took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
                    // printf("dotVector\n");
                    // gettimeofday(&start, NULL);
                    double res = m->bridges[index].dotVector_f(row->extra, col->extra, chunk_size);
                    // gettimeofday(&stop, NULL);
                    // printf("   dot product took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
                    // printf("%d %d %d res=%f\n", j, c, i, res);
                    r->value.f[i * r->ncol + j] += res;
                    m->bridges[index].vector_delete(row);
                    free(A_row);
                }
                m->bridges[index].vector_delete(col);
                free(B_col);
            }
            printf("j=%d ", j);
            gettimeofday(&end, NULL);
            timersub(&end, &ini, &tval_result);
            printf("Time elapsed: %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
            // printf("col took %lu us\n\n", (end.tv_sec - ini.tv_sec) * 1000000 + end.tv_usec - ini.tv_usec);
            // getchar();
        }
    }

    // clear_input(i, 2);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

void **vec_mulsc(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    double scalar = 1.0;
    if (type(*in[0]) == T_FLOAT)
    {
        scalar = fvalue(*in[0]);
    }
    else if (type(*in[0]) == T_INT)
    {
        scalar = (double)ivalue(*in[0]);
    }
    else
    {
        runerror("Runtime error: invalid type\n");
    }
    vector_t *v = (vector_t *)vvalue(*in[1]);
    m->bridges[index].vecreqdev(v);

    vector_t *r = m->bridges[index].vector_new(v->len, T_FLOAT, 0, NULL);
    r->location = LOCDEV;

    r->extra = (void *)m->bridges[index].mulScalarVector_f(v->extra, scalar, v->len);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

vector_t *vec_mul_complex_scalar(bridge_manager_t *m, int index, complex_t *s, vector_t *a)
{

    m->bridges[index].vecreqdev(a);

    vector_t *r = m->bridges[index].vector_new(a->len, T_COMPLEX, 0, NULL); //(vector_t *) malloc( sizeof( vector_t ) );
    r->location = LOCDEV;
    // apenas cpu
    //  free( r->value.f);
    r->extra = (void *)m->bridges[index].mulComplexScalarVector_f(a->extra, s->re, s->im, a->len);

    return (void *)r;
}

vector_t *mul_complex_scalar_complex_vec(bridge_manager_t *m, int index, complex_t *s, vector_t *a)
{
    m->bridges[index].vecreqdev(a);

    vector_t *r = m->bridges[index].vector_new(a->len, T_COMPLEX, 0, NULL); //(vector_t *) malloc( sizeof( vector_t ) );
    r->location = LOCDEV;
    // apenas cpu
    //  free( r->value.f);
    r->extra = (void *)m->bridges[index].mulComplexScalarComplexVector_f(a->extra, s->re, s->im, a->len);

    return (void *)r;
}

vector_t *mul_float_scalar_complex_vec(bridge_manager_t *m, int index, double d, vector_t *a)
{
    m->bridges[index].vecreqdev(a);

    vector_t *r = m->bridges[index].vector_new(a->len, T_COMPLEX, 0, NULL); //(vector_t *) malloc( sizeof( vector_t ) );
    r->location = LOCDEV;
    // apenas cpu
    //  free( r->value.f);

    r->extra = (void *)m->bridges[index].mulFloatScalarComplexVector_f(a->extra, d, a->len);

    return (void *)r;
}

void **mat_mulsc(bridge_manager_t *mg, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    double scalar = 1.0;
    if (type(*in[0]) == T_FLOAT)
    {
        scalar = fvalue(*in[0]);
    }
    else if (type(*in[0]) == T_INT)
    {
        scalar = (double)ivalue(*in[0]);
    }
    else
    {
        runerror("Runtime error: invalid type\n");
    }

    matrix_t *m = (matrix_t *)vvalue(*in[1]);
    mg->bridges[index].matreqdev(m);

    matrix_t *r = NULL;
    if (m->type == T_FLOAT)
    {
        r = mg->bridges[index].matrix_new(m->nrow, m->ncol, T_FLOAT, 0, NULL);
        r->extra = mg->bridges[index].mulScalarVector_f(m->extra, scalar, m->nrow * m->ncol);
        r->location = LOCDEV;
    }
    else if (m->type == T_COMPLEX)
    {
        r = mg->bridges[index].matrix_new(m->nrow, m->ncol, T_COMPLEX, 0, NULL);
        r->extra = mg->bridges[index].mulScalarVector_f(m->extra, scalar, 2 * m->nrow * m->ncol);
        r->location = LOCDEV;
    }

    clear_input(i, 2);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

matrix_t *mul_complex_scalar_complex_mat(bridge_manager_t *mg, int index, complex_t *s, matrix_t *m)
{
    matrix_t *r = NULL;
    mg->bridges[index].matreqdev(m);
    r = mg->bridges[index].matrix_new(m->nrow, m->ncol, T_COMPLEX, 0, NULL);
    r->extra = mg->bridges[index].mulComplexScalarComplexVector_f(m->extra, s->re, s->im, m->nrow * m->ncol);
    r->location = LOCDEV;

    return (void *)r;
}

matrix_t *mul_complex_scalar_float_mat(bridge_manager_t *mg, int index, complex_t *s, matrix_t *m)
{
    matrix_t *r = NULL;
    mg->bridges[index].matreqdev(m);
    r = mg->bridges[index].matrix_new(m->nrow, m->ncol, T_COMPLEX, 0, NULL);
    r->extra = mg->bridges[index].mulComplexScalarVector_f(m->extra, s->re, s->im, m->nrow * m->ncol);
    r->location = LOCDEV;

    return (void *)r;
}

void **matvec_mul3(bridge_manager_t *mg, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    vector_t *v = (vector_t *)vvalue(*in[0]);
    // vector_t * r = (vector_t *) malloc( sizeof( vector_t ) );
    vector_t *r;

    // do I have to assume that it needs to be copied everytime?
    if (v->location != LOCDEV)
    {
        mg->bridges[index].vecreqdev(v);
    }

    if (type(*in[1]) == T_MATRIX)
    {
        matrix_t *m = (matrix_t *)vvalue(*in[1]);
        r = mg->bridges[index].vector_new(m->nrow, m->type, 0, NULL);
        mg->bridges[index].vecreqdev(r);
        if (m->location != LOCDEV)
        {
            mg->bridges[index].matreqdev(m);
        }

        if (m->type == T_FLOAT && v->type == T_FLOAT)
        {
            r->extra = (void *)mg->bridges[index].matVecMul3_f(m->extra, v->extra, m->ncol, m->nrow);
        }
        else if (m->type == T_COMPLEX && v->type == T_COMPLEX)
        {
            r->extra = (void *)mg->bridges[index].matVecMul3Complex_f(m->extra, v->extra, m->ncol, m->nrow);
        }
        if (status != NULL)
        {
            *status = 0;
        }
        return (void *)r;
    }
    else if (type(*in[1]) == T_SMATRIX)
    {
        smatrix_t *m = (smatrix_t *)vvalue(*in[1]);
        r = mg->bridges[index].vector_new(m->nrow, m->type, 0, NULL);
        mg->bridges[index].vecreqdev(r);
        if (m->location != LOCDEV)
        {
            mg->bridges[index].smatreqdev(m);
        }
        if (m->type == T_FLOAT && v->type == T_FLOAT)
        {
            r->extra = (void *)mg->bridges[index].sparseVecMul_f(m->extra, m->idxColMem, v->extra, m->nrow, m->maxcols);
        }
        else if (m->type == T_COMPLEX && v->type == T_COMPLEX)
        {
            r->extra = (void *)mg->bridges[index].sparseComplexVecMul_f(m->extra, m->idxColMem, v->extra, m->nrow, m->maxcols);
        }
        if (status != NULL)
        {
            *status = 0;
        }
        return (void *)r;
    }
    else if (type(*in[1]) == T_RMATRIX)
    {
        if (status != NULL)
        {
            *status = -1;
        }
        return (void **)NULL;
    }
    else
    {
        if (status != NULL)
        {
            *status = -1;
        }
        return (void **)NULL;
    }
}

object_t **convert_int_and_vector_to_object(int int_value, vector_t *vector_a)
{
    object_t **in;
    in = (object_t **)malloc(2 * sizeof(object_t *));

    in[0] = (object_t *)malloc(sizeof(object_t));
    ivalue(*in[0]) = int_value;
    in[0]->type = T_INT;

    in[1] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[1]) = vector_a;
    in[1]->type = T_VECTOR;

    return in;
}

void **vec_add_off(bridge_manager_t *m, int index, void **i, int *status)
{
    object_t **in = (object_t **)i;
    int offset = ivalue(*in[0]);
    vector_t *a = (vector_t *)vvalue(*in[1]);
    int parts = a->len / offset;

    vector_t *r = m->bridges[index].vector_new(offset, T_FLOAT, 0, NULL);
    // apenas cpu
    //  free( r->value.f);
    m->bridges[index].vecreqdev(a);
    m->bridges[index].vecreqdev(r);

    r->extra = (void *)m->bridges[index].vecAddOff_f(a->extra, offset, parts);

    clear_input(i, 2);
    if (status != NULL)
    {
        *status = 0;
    }
    return (void *)r;
}

// Converts

object_t **convert_vectors_to_object(vector_t *vector_a, vector_t *vector_b)
{
    object_t **in;
    if (vector_b != NULL)
    {
        in = (object_t **)malloc(2 * sizeof(object_t *));
        in[1] = (object_t *)malloc(sizeof(object_t));
        in[1]->value.v = (void *)vector_b;
        in[1]->type = T_VECTOR;
    }
    else
    {
        in = (object_t **)malloc(2 * sizeof(object_t *));
    }

    in[0] = (object_t *)malloc(sizeof(object_t));
    in[0]->value.v = (void *)vector_a;
    in[0]->type = T_VECTOR;
    return in;
}

object_t **convert_vector_and_matrix_to_object(vector_t *vector_a, matrix_t *vector_b)
{
    object_t **in;
    if (vector_b != NULL)
    {
        in = (object_t **)malloc(2 * sizeof(object_t *));
        in[1] = (object_t *)malloc(sizeof(object_t));
        vvalue(*in[1]) = vector_b;
        in[1]->type = T_MATRIX;
    }
    else
    {
        in = (object_t **)malloc(sizeof(object_t *));
    }

    in[0] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[0]) = vector_a;
    in[0]->type = T_VECTOR;

    return in;
}

object_t **convert_matrices_to_object(matrix_t *matrix_a, matrix_t *matrix_b)
{
    object_t **in;
    in = (object_t **)malloc(2 * sizeof(object_t *));

    in[0] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[0]) = matrix_a;
    in[0]->type = T_MATRIX;

    in[1] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[1]) = matrix_b;
    in[1]->type = T_MATRIX;

    return in;
}

object_t **convert_scalar_and_vector_to_object(double s, vector_t *vector_a)
{
    object_t **in;
    in = (object_t **)malloc(2 * sizeof(object_t *));

    in[0] = (object_t *)malloc(sizeof(object_t));
    fvalue(*in[0]) = s;
    in[0]->type = T_FLOAT;

    in[1] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[1]) = vector_a;
    in[1]->type = T_VECTOR;

    return in;
}

object_t **convert_scalar_and_matrix_to_object(double scalar, matrix_t *matrix_a)
{
    object_t **in;
    in = (object_t **)malloc(2 * sizeof(object_t *));

    in[0] = (object_t *)malloc(sizeof(object_t));
    fvalue(*in[0]) = scalar;
    in[0]->type = T_FLOAT;

    in[1] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[1]) = matrix_a;
    in[1]->type = T_MATRIX;

    return in;
}

object_t **convert_vector_and_smatrix_to_object(vector_t *vector_a, smatrix_t *smatrix_b)
{
    object_t **in;
    if (smatrix_b != NULL)
    {
        in = (object_t **)malloc(2 * sizeof(object_t *));
        in[1] = (object_t *)malloc(sizeof(object_t));
        vvalue(*in[1]) = smatrix_b;
        in[1]->type = T_SMATRIX;
    }
    else
    {
        in = (object_t **)malloc(sizeof(object_t *));
    }

    in[0] = (object_t *)malloc(sizeof(object_t));
    vvalue(*in[0]) = vector_a;
    in[0]->type = T_VECTOR;

    return in;
}
