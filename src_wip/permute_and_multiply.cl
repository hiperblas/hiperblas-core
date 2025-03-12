__kernel void permute_and_multiply(
    __global const int* A_col_idx,  // Column indices (mapping rows in A)
    __global const int* B_row_ptr,  // Row pointer of block diagonal matrix B
    __global const int* B_col_idx,  // Column indices of B
    __global const float* B_values, // Values of B
    __global int* C_row_ptr,        // Row pointer for C
    __global int* C_col_idx,        // Column indices for C
    __global float* C_values,       // Values for C
    int A_nrows
) {
    int row = get_global_id(0); // Each thread processes one row of A
    
    int thread = 0;
    if(row == thread)
    {
        printf("Thread %d:\n", row);
        printf("%d \n", A_nrows);
    }
    
    if (row < A_nrows) {
        // Find the column index (row index in B) where A has its nonzero
        int permuted_row = A_col_idx[row]; // Since A is a permutation matrix

        // Copy row `permuted_row` from B into row `row` in C
        int startB = B_row_ptr[permuted_row];
        int endB = B_row_ptr[permuted_row + 1];


        int startC = C_row_ptr[row]; // Start index for C

        for (int i = 0; i < (endB - startB); i++) {
            C_col_idx[startC + i] = B_col_idx[startB + i];  // Copy column index
            C_values[startC + i] = B_values[startB + i];    // Copy value
            if(row == thread)
            {
                printf("%d: %f\n", C_col_idx[startC + i], C_values[startC + i]);
            }
        }
    }
}