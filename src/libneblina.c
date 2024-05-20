#include "libneblina.h"
#include <stdio.h>
#include <stdlib.h>
#include <bridge_api.h>
#include <dlfcn.h>

void load_function(bridge_manager_t *manager, void *(**function_ptr)(), char *function_name, int index)
{
    void *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_double_function(bridge_manager_t *manager, double (**function_ptr)(), char *function_name, int index)
{
    double (*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_int_function(bridge_manager_t *manager, int (**function_ptr)(), char *function_name, int index)
{
    int (*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_void_function(bridge_manager_t *manager, void (**function_ptr)(), char *function_name, int index)
{
    void (*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_vector_function(bridge_manager_t *manager, vector_t *(**function_ptr)(), char *function_name, int index)
{
    vector_t *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_matrix_function(bridge_manager_t *manager, matrix_t *(**function_ptr)(), char *function_name, int index)
{
    matrix_t *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_smatrix_function(bridge_manager_t *manager, smatrix_t *(**function_ptr)(), char *function_name, int index)
{
    smatrix_t *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_slist_function(bridge_manager_t *manager, slist *(**function_ptr)(), char *function_name, int index)
{
    slist *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_complex_function(bridge_manager_t *manager, complex_t *(**function_ptr)(), char *function_name, int index)
{
    complex_t *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_double_pointer_function(bridge_manager_t *manager, double *(**function_ptr)(), char *function_name, int index)
{
    double *(*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

void load_long_function(bridge_manager_t *manager, long (**function_ptr)(), char *function_name, int index)
{
    long (*externalFunction)(void);
    *(void **)(&externalFunction) = dlsym(manager->bridges[index].plugin_handle, function_name);
    *function_ptr = externalFunction;
    char *result = dlerror();
    if (result)
    {
        printf("Cannot find init in %s: %s", function_name, result);
    }
}

/*
*
 TO DO
 * como passar a pasta de instalacao dos plugins
 * como passar o nome dos plugins existentes
 * como passar do python para o C o indice do dispositivo
 * como fazer testes unitarios com as funcoes de cada plugin carregado
 * * no opencl bridge esta ok ja, precisa criar uma implementacao do cpu bridge
 * como implementar o release plugins
 *
 */

void load_plugin(bridge_manager_t *manager, char *library_name, int index)
{
    char *plugin_name = NULL;

    manager->bridges[index].plugin_handle = dlopen(library_name, RTLD_NOW);
    if (!manager->bridges[index].plugin_handle)
    {
        if (plugin_name != NULL)
        {
            printf("Cannot load %s: %s", plugin_name, dlerror());
        }
    }

    load_function(manager, &(manager->bridges[index].addVectorF_f), "addVectorF", index);

    load_function(manager, &(manager->bridges[index].addVectorC_f), "addVectorC", index);

    load_function(manager, &(manager->bridges[index].addVectorFC_f), "addVectorFC", index);

    load_function(manager, &(manager->bridges[index].vecAddOff_f), "vecAddOff", index);

    load_function(manager, &(manager->bridges[index].mulScalarVector_f), "mulScalarVector", index);

    load_function(manager, &(manager->bridges[index].subVector_f), "subVector", index);

    load_function(manager, &(manager->bridges[index].subVectorC_f), "subVectorC", index);

    load_function(manager, &(manager->bridges[index].vecConjugate_f), "vecConjugate", index);

    load_function(manager, &(manager->bridges[index].prodVector_f), "prodVector", index);

    load_function(manager, &(manager->bridges[index].prodComplexVector_f), "prodComplexVector", index);

    load_double_function(manager, &(manager->bridges[index].sumVector_f), "sumVector", index);

    load_double_function(manager, &(manager->bridges[index].normVector_f), "normVector", index);

    load_double_function(manager, &(manager->bridges[index].dotVector_f), "dotVector", index);

    load_void_function(manager, &(manager->bridges[index].dotVectorComplex_f), "dotVectorComplex", index);

    load_vector_function(manager, &(manager->bridges[index].vector_new), "vector_new", index);

    load_void_function(manager, &(manager->bridges[index].vector_delete), "vector_delete", index);

    load_void_function(manager, &(manager->bridges[index].vecreqdev), "vecreqdev", index);

    load_void_function(manager, &(manager->bridges[index].vecreqhost), "vecreqhost", index);

    load_void_function(manager, &(manager->bridges[index].InitEngine_f), "InitEngine", index);

    load_void_function(manager, &(manager->bridges[index].StopEngine_f), "StopEngine", index);

    load_long_function(manager, &(manager->bridges[index].get_Engine_Max_Memory_Allocation_f), "get_Engine_Max_Memory_Allocation", index);

    load_int_function(manager, &(manager->bridges[index].list_len), "list_len", index);

    load_complex_function(manager, &(manager->bridges[index].complex_new), "complex_new", index);

    load_void_function(manager, &(manager->bridges[index].complex_delete), "complex_delete", index);

    load_function(manager, &(manager->bridges[index].mulComplexScalarVector_f), "mulComplexScalarVector", index);

    load_function(manager, &(manager->bridges[index].mulComplexScalarComplexVector_f), "mulComplexScalarComplexVector", index);

    load_function(manager, &(manager->bridges[index].mulFloatScalarComplexVector_f), "mulFloatScalarComplexVector", index);

    load_matrix_function(manager, &(manager->bridges[index].matrix_new), "matrix_new", index);

    load_void_function(manager, &(manager->bridges[index].matrix_delete), "matrix_delete", index);

    load_double_function(manager, &(manager->bridges[index].matrix_get_complex_imaginary_value), "matrix_get_complex_imaginary_value", index);

    load_double_function(manager, &(manager->bridges[index].matrix_get_complex_real_value), "matrix_get_complex_real_value", index);

    load_double_function(manager, &(manager->bridges[index].matrix_get_real_value), "matrix_get_real_value", index);

    load_void_function(manager, &(manager->bridges[index].matrix_set_complex_value), "matrix_set_complex_value", index);

    load_double_pointer_function(manager, &(manager->bridges[index].matrix_copy_col), "matrix_copy_col", index);

    load_double_pointer_function(manager, &(manager->bridges[index].matrix_copy_row), "matrix_copy_row", index);

    load_void_function(manager, &(manager->bridges[index].matrix_set_real_value), "matrix_set_real_value", index);

    load_void_function(manager, &(manager->bridges[index].matreqdev), "matreqdev", index);

    load_void_function(manager, &(manager->bridges[index].matreqhost), "matreqhost", index);

    load_function(manager, &(manager->bridges[index].matMul_f), "matMul", index);

    load_function(manager, &(manager->bridges[index].matVecMul3_f), "matVecMul3", index);

    load_function(manager, &(manager->bridges[index].matVecMul3Complex_f), "matVecMul3Complex", index);

    load_smatrix_function(manager, &(manager->bridges[index].smatrix_new), "smatrix_new", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_t_clear), "smatrix_t_clear", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_load_double), "smatrix_load_double", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_set_real_value), "smatrix_set_real_value", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_pack), "smatrix_pack", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_set_complex_value), "smatrix_set_complex_value", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_pack_complex), "smatrix_pack_complex", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_load_complex), "smatrix_load_complex", index);

    load_void_function(manager, &(manager->bridges[index].smatrix_delete), "smatrix_delete", index);

    load_void_function(manager, &(manager->bridges[index].smatreqhost), "smatreqhost", index);

    load_void_function(manager, &(manager->bridges[index].smatreqdev), "smatreqdev", index);

    load_slist_function(manager, &(manager->bridges[index].slist_add), "slist_add", index);

    load_void_function(manager, &(manager->bridges[index].slist_clear), "slist_clear", index);

    load_function(manager, &(manager->bridges[index].sparseVecMul_f), "sparseVecMul", index);

    load_function(manager, &(manager->bridges[index].sparseComplexVecMul_f), "sparseComplexVecMul", index);
}

void release_plugin(bridge_manager_t *manager, int index)
{
    // TODO currently we cannot close the library because the Python interpreter still needs the
    // deletion functions

    //    int error = dlclose(manager->bridges[index].plugin_handle);
    //
    //    if (error) {
    //        printf("Cannot load : %s", dlerror());
    //    }
}

void ord_smat(double *m, int *index, int max, int N)
{
    int i, j, tmpi, k;
    double tmpd;
    for (k = 0; k < N; k++)
    {
        for (i = max - 1; i >= 1; i--)
        {
            for (j = 0; j < i; j++)
            {
                if (index[k * max + j + 1] != -1 && index[k * max + j] != -1 && index[k * max + j] > index[k * max + j + 1])
                {
                    tmpi = index[k * max + j];
                    index[k * max + j] = index[k * max + j + 1];
                    index[k * max + j + 1] = tmpi;

                    tmpd = m[2 * (k * max + j)];
                    m[2 * (k * max + j)] = m[2 * (k * max + j + 1)];
                    m[2 * (k * max + j + 1)] = tmpd;

                    tmpd = m[2 * (k * max + j) + 1];
                    m[2 * (k * max + j) + 1] = m[2 * (k * max + j + 1) + 1];
                    m[2 * (k * max + j + 1) + 1] = tmpd;
                }
            }
        }
    }
}

void smatrix_line_to_col(double *out, int *index_out, double *in, int *index_in, int max, int N)
{
    int i, col, j;
    int count[N];
    for (i = 0; i < N; i++)
    {
        count[i] = 0;
        for (j = 0; j < max; j++)
            index_out[i * max + j] = -1;
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < max; j++)
        {
            col = index_in[i * max + j];
            index_out[col * max + count[col]] = i;
            out[2 * (col * max + count[col])] = in[2 * (i * max + j)];
            out[2 * (col * max + count[col]) + 1] = in[2 * (i * max + j) + 1];
            count[col]++;
        }
    }
}

void print_data_type(data_type type)
{
    switch (type)
    {
    case T_STRING:
        printf("type: T_STRING\n");
        break;
    case T_INT:
        printf("type: T_INT\n");
        break;
    case T_FLOAT:
        printf("type: T_FLOAT\n");
        break;
    case T_ADDR:
        printf("type: T_ADDR\n");
        break;
    case T_NDEF:
        printf("type: T_NDEF\n");
        break;
    case T_LIST:
        printf("type: T_LIST\n");
        break;
    case T_STRTOREL:
        printf("type: T_STRTOREL\n");
        break;
    case T_CFUNC:
        printf("type: T_CFUNC\n");
        break;
    case T_VECTOR:
        printf("type: T_VECTOR\n");
        break;
    case T_MATRIX:
        printf("type: T_MATRIX\n");
        break;
    case T_SMATRIX:
        printf("type: T_SMATRIX\n");
        break;
    case T_COMPLEX:
        printf("type: T_COMPLEX\n");
        break;
    case T_FILE:
        printf("type: T_FILE\n");
        break;
    case T_ANY:
        printf("type: T_ANY\n");
        break;
    default:
        printf("type: INDEFINED\n");
    };
}

void clear_input(void **input, int n_params)
{
    int k = 0;
    object_t **in = (object_t **)input;
    for (k = 0; k < n_params; k++)
    {
        if (type(*in[k]) == T_STRING)
            free(svalue(*in[k]));

        if (type(*in[k]) == T_COMPLEX)
        {
            complex_t *r = (complex_t *)vvalue(*in[k]);
            free(r);
        }
    }
}

void neblina_strtype(data_type type, char out[256])
{
    switch (type)
    {
    case (T_STRING):
        sprintf(out, "string");
        break;
    case (T_INT):
        sprintf(out, "int");
        break;
    case (T_FLOAT):
        sprintf(out, "double");
        break;
    case (T_VECTOR):
        sprintf(out, "vector");
        break;
    case (T_MATRIX):
        sprintf(out, "matrix");
        break;
    case (T_FILE):
        sprintf(out, "file");
        break;
    case (T_COMPLEX):
        sprintf(out, "complex");
        break;
    case (T_LIST):
        sprintf(out, "list");
        break;
    default:
        sprintf(out, "unknown:%d", type);
        break;
    }
}
