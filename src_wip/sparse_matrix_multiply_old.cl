inline void atomic_add_float(volatile __local float *address, float value) {
    union {
        unsigned int int_value;
        float float_value;
    } old, new_val;

    do {
        old.float_value = *address;
        new_val.float_value = old.float_value + value;
    } while (atomic_cmpxchg((volatile __local unsigned int *)address,
                            old.int_value,
                            new_val.int_value) != old.int_value);
}

__kernel void sparse_matrix_multiply_csr_column_major(
    __global const int *row_ptr_a,
    __global const int *col_idx_a, 
    __global const float *values_a,
    __global const int *col_ptr_b,  // Column pointers (column-major CSR for B)
    __global const int *row_idx_b, // Row indices for column-major CSR
    __global const float *values_b, // Non-zero values for B in column-major
    __global int *row_ptr_c, 
    __global int *col_idx_c, 
    __global float *values_c,
    const int nrows_a,
    const int ncols_b) {
    
    int row = get_global_id(0); // Each work-item processes one row of A
    if (row >= nrows_a) return;

    // Allocate local memory to temporarily store values for row C
    __local float temp_values[1024]; // Adjustable size
    __local int temp_cols[1024];     // Adjustable size
    __local int temp_count;

    if (get_local_id(0) == 0) {
        temp_count = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Iterate over the elements of the current row of A
    for (int i = row_ptr_a[row]; i < row_ptr_a[row + 1]; i++) {
        int col_a = col_idx_a[i];    // Column index from A
        float val_a = values_a[i];  // Value from A

        // Traverse the corresponding column of B (column-major CSR)
        for (int j = col_ptr_b[col_a]; j < col_ptr_b[col_a + 1]; j++) {
            int row_b = row_idx_b[j]; // Row index from column-major B
            float val_b = values_b[j]; // Value from B

            // Search for the corresponding index in temp_cols to accumulate the result
            int found = 0;
            for (int k = 0; k < temp_count; k++) {
                if (temp_cols[k] == row_b) {
                    atomic_add_float(&temp_values[k], val_a * val_b);
                    found = 1;
                    break;
                }
            }

            // If not found, add a new entry
            if (!found) {
                int idx = atomic_add(&temp_count, 1);
                temp_cols[idx] = row_b;
                temp_values[idx] = val_a * val_b;
            }
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // Write the accumulated values back to C
    int start_idx = atomic_add(&row_ptr_c[row], temp_count);
    for (int k = 0; k < temp_count; k++) {
        col_idx_c[start_idx + k] = temp_cols[k];
        values_c[start_idx + k] = temp_values[k];
    }
}