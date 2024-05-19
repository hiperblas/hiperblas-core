#ifndef __NEBLINASTD___
#define __NEBLINASTD___
#ifdef	__cplusplus
extern "C" {
#endif

#include "neblina.h"
#include "neblina_vector.h"
#include "neblina_matrix.h"
#include "neblina_smatrix.h"
#include "neblina_complex.h"
#include "bridge_api.h"
    
#define new_str(i) (char *) malloc( sizeof(char)*(i))

 void delete_object_array(object_t ** in, int len);
 void delete_object(object_t * in);
 
 object_t ** convert_vectors_to_object(vector_t * vector_a, vector_t * vector_b);
 object_t ** convert_matrices_to_object(matrix_t * matrix_a, matrix_t * matrix_b);
 object_t **convert_int_and_vector_to_object(int int_value, vector_t *vector_a);
 object_t ** convert_vector_and_matrix_to_object(vector_t * vector_a, matrix_t * matrix_b);
 object_t ** convert_vector_and_smatrix_to_object(vector_t * vector_a, smatrix_t * smatrix_b);
 object_t ** convert_scalar_and_vector_to_object(double scalar, vector_t * vector_a);
 object_t ** convert_scalar_and_matrix_to_object(double scalar, matrix_t * matrix_a);

 void ** neblina_type   ( void ** i, int * status );
 int     vec_len        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_add        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_prod       ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_conj       ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_sub        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_add_off    ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_add_off2   ( void ** i, int * status );
 void ** vec_mulsc      ( bridge_manager_t *m, int index, void ** i, int * status ); 
 vector_t * vec_mul_complex_scalar ( bridge_manager_t *m, int index, complex_t * s, vector_t * a); 
 vector_t * mul_complex_scalar_complex_vec( bridge_manager_t *m, int index, complex_t * s, vector_t * a); 
 vector_t * mul_float_scalar_complex_vec( bridge_manager_t *m, int index, double d, vector_t * a); 
 void ** vec_sum        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_norm       ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** vec_dot        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** mat_len_col    ( void ** i, int * status );
 void ** mat_len_row    ( void ** i, int * status );
 void ** mat_mul        ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** mat_add        ( bridge_manager_t *m, int index, void ** i, int * status ); 
 void ** mat_sub        ( bridge_manager_t *m, int index, void ** i, int * status );
 matrix_t * mul_complex_scalar_complex_mat( bridge_manager_t *m, int index, complex_t * s, matrix_t * a);
 matrix_t * mul_complex_scalar_float_mat( bridge_manager_t *m, int index, complex_t * s, matrix_t * a);
 void ** mat_mulsc      ( bridge_manager_t *m, int index, void ** i, int * status ); 
 void ** mat_mulscrow   ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** mat_mulsccol   ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** mat_transp     ( bridge_manager_t *m, int index, void ** i, int * status );
 void ** matvec_mul1    ( void ** i, int * status );
 void ** matvec_mul2    ( void ** i, int * status );
 void ** matvec_mul3    ( bridge_manager_t *m, int index, void ** i, int * status );

 void ** vec_sub_cpu    ( void ** i, int * status );
 void ** vec_add_cpu    ( void ** i, int * status );
 void ** mat_mul_cpu    ( void ** i, int * status );
 void ** mat_add_cpu    ( void ** i, int * status );
 void ** mat_sub_cpu    ( void ** i, int * status );
 void ** vec_mulsc_cpu  ( void ** i, int * status );
 void ** mat_mulsc_cpu  ( void ** i, int * status );
 void ** matvec_mul_cpu ( void ** i, int * status );
 void ** smatvec_mul_cpu( void ** i, int * status );
 void ** mat_transp_cpu ( void ** i, int * status );
 void ** vec_sum_cpu    ( void ** i, int * status );
 void ** vec_norm_cpu   ( void ** i, int * status );
 void ** vec_dot_cpu    ( void ** i, int * status );
 void ** toint          ( void ** i, int * status );
 void ** tostr          ( void ** i, int * status );
 void ** tostr2         ( void ** i, int * status );
 void ** neblina_at     ( void ** i, int * status );
 void ** neblina_upper  ( void ** i, int * status );
 void ** neblina_lower  ( void ** i, int * status );
 void ** todouble       ( void ** i, int * status );
 //void ** complex_new    ( void ** i, int * status );
 void ** complex_real   ( void ** i, int * status );
 void ** complex_imag   ( void ** i, int * status );
 void ** complex_conj   ( void ** i, int * status );
 void ** neblina_ludecomp( void ** i, int * status );
 void ** neblina_solve( void ** i, int * status );
 void ** neblina_list_new( void ** i, int * status );
 void ** neblina_list_append( void ** i, int * status );
 void ** neblina_list_get( void ** i, int * status );


 void ** neblina_mat_sqr( void ** i, int * s );

void ** init     ( void ** i, int * status );

 void ** neblina_rmatrix( void ** i, int * status );

  #ifdef	__cplusplus
}
#endif

#endif
