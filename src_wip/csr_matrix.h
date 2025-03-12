#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __csr_t {
    int nrow;
    int ncol;
    int nnz;

    int* row_ptr;
    int* col_idx;
    float* values;
} csr_t;

#ifdef __cplusplus
}
#endif

#endif // CSR_MATRIX_H