#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>

// OpenCL Kernel for CSR matrix multiplication
const char* csr_matrix_mult_kernel = R"CLC(
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
)CLC";

// Function to check OpenCL errors
void check_error(cl_int err, const char* operation) {
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Error during %s: %d\n", operation, err);
        exit(EXIT_FAILURE);
    }
}

// Host function to compute C_row_ptr
void compute_C_row_ptr(const int* A_col_idx, const int* B_row_ptr, int* C_row_ptr, int A_nrows) {
    C_row_ptr[0] = 0; // First row starts at index 0
    for (int i = 0; i < A_nrows; i++) {
        int permuted_row = A_col_idx[i];  // Get the row index in B
        int nnz_in_B_row = B_row_ptr[permuted_row + 1] - B_row_ptr[permuted_row];
        C_row_ptr[i + 1] = C_row_ptr[i] + nnz_in_B_row;
    }
}

// Main function
int main() {
    // Matrices A and B in CSR format
    /*
        [ 0  1  0  0 ]
        [ 1  0  0  0 ]
        [ 0  0  0  1 ]
        [ 0  0  1  0 ]
    */
    int A_row_ptr[] = {0, 1, 2, 3, 4};
    int A_col_idx[] = {1, 0, 3, 2};
    float A_values[] = {1.0, 1.0, 1.0, 1.0};
    int A_nrows = 4;
    int A_ncols = 4;

    /*
        [ 1  2  0  0 ]
        [ 2  1  0  0 ]
        [ 0  0  3  4 ]
        [ 0  0  4  3 ]
    */
    int B_row_ptr[] = {0, 2, 4, 6, 8};
    int B_col_idx[] = {0, 1, 0, 1, 2, 3, 2, 3};
    float B_values[] = {1.0, 2.0, 2.0, 1.0, 3.0, 4.0, 4.0, 3.0};
    int B_nrows = 4;
    int B_ncols = 4;

    /*
        [ 2  1  0  0 ]
        [ 1  2  0  0 ]
        [ 0  0  4  3 ]
        [ 0  0  3  4 ]
    */
    // Result matrix C (dynamic allocation)
    int C_row_ptr[A_nrows + 1];
    int C_col_idx[1024]; // Large allocation for simplicity
    float C_values[1024];

    // Initialize C_row_ptr
    compute_C_row_ptr(A_col_idx, B_row_ptr, C_row_ptr, A_nrows);
    /*
    for (int i = 0; i <= A_nrows; i++) {
        printf("%d ", B_row_ptr[i]);
    }
    */
    
    // OpenCL initialization
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;
    cl_int err;

    // Get platform and device
    err = clGetPlatformIDs(1, &platform, NULL);
    check_error(err, "clGetPlatformIDs");

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    check_error(err, "clGetDeviceIDs");

    // Create context and command queue
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    check_error(err, "clCreateContext");

    queue = clCreateCommandQueue(context, device, 0, &err);
    check_error(err, "clCreateCommandQueue");

    // Create program and kernel
    program = clCreateProgramWithSource(context, 1, &csr_matrix_mult_kernel, NULL, &err);
    check_error(err, "clCreateProgramWithSource");

    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        char buffer[2048];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
        fprintf(stderr, "Error building program:\n%s\n", buffer);
        exit(EXIT_FAILURE);
    }

    kernel = clCreateKernel(program, "permute_and_multiply", &err);
    check_error(err, "clCreateKernel");

    // Create buffers
    //cl_mem A_row_ptr_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(A_row_ptr), A_row_ptr, &err);
    //check_error(err, "clCreateBuffer A_row_ptr");

    cl_mem A_col_idx_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(A_col_idx), A_col_idx, &err);
    check_error(err, "clCreateBuffer A_col_idx");

    //cl_mem A_values_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(A_values), A_values, &err);
    //check_error(err, "clCreateBuffer A_values");

    cl_mem B_row_ptr_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(B_row_ptr), B_row_ptr, &err);
    check_error(err, "clCreateBuffer B_row_ptr");

    cl_mem B_col_idx_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(B_col_idx), B_col_idx, &err);
    check_error(err, "clCreateBuffer B_col_idx");

    cl_mem B_values_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(B_values), B_values, &err);
    check_error(err, "clCreateBuffer B_values");

    //cl_mem C_row_ptr_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(C_row_ptr), NULL, &err);
    //check_error(err, "clCreateBuffer C_row_ptr");
    cl_mem C_row_ptr_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(C_row_ptr), C_row_ptr, &err);
    check_error(err, "clCreateBuffer C_row_ptr");

    cl_mem C_col_idx_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(C_col_idx), NULL, &err);
    check_error(err, "clCreateBuffer C_col_idx");

    cl_mem C_values_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(C_values), NULL, &err);
    check_error(err, "clCreateBuffer C_values");

    // Set kernel arguments
    //__global const int* A_row_ptr,
    //err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &A_row_ptr_buf);
    //check_error(err, "clSetKernelArg 0");

    //__global const int* A_col_idx,
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &A_col_idx_buf);
    check_error(err, "clSetKernelArg 0");

    //__global const float* A_values,
    //err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &A_values_buf);
    //check_error(err, "clSetKernelArg 2");

    //__global const int* B_row_ptr,
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &B_row_ptr_buf);
    check_error(err, "clSetKernelArg 1");

    //__global const int* B_col_idx,
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &B_col_idx_buf);
    check_error(err, "clSetKernelArg 2");

    //__global const float* B_values,
    err = clSetKernelArg(kernel, 3, sizeof(cl_mem), &B_values_buf);
    check_error(err, "clSetKernelArg 3");

    //__global int* C_row_ptr,
    err = clSetKernelArg(kernel, 4, sizeof(cl_mem), &C_row_ptr_buf);
    check_error(err, "clSetKernelArg 4");

    //__global int* C_col_idx,
    err = clSetKernelArg(kernel, 5, sizeof(cl_mem), &C_col_idx_buf);
    check_error(err, "clSetKernelArg 5");

    //__global float* C_values,
    err = clSetKernelArg(kernel, 6, sizeof(cl_mem), &C_values_buf);
    check_error(err, "clSetKernelArg 6");

    //printf("> %d %d", A_nrows, B_nrows); //OK

    //const int A_nrows,
    err = clSetKernelArg(kernel, 7, sizeof(int), &A_nrows);
    check_error(err, "clSetKernelArg 7");

    //const int A_ncols
    //err = clSetKernelArg(kernel, 10, sizeof(int), &A_ncols);
    //check_error(err, "clSetKernelArg 10");

    //const int B_ncols
    //err = clSetKernelArg(kernel, 11, sizeof(int), &B_ncols);
    //check_error(err, "clSetKernelArg 11");

    size_t global_work_size = A_nrows;
    //printf("> %d", global_work_size);
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
    check_error(err, "clEnqueueNDRangeKernel");

    // Wait for the kernel to complete
    clFinish(queue);

    // Read back the results
    //err = clEnqueueReadBuffer(queue, C_row_ptr_buf, CL_TRUE, 0, sizeof(C_row_ptr), C_row_ptr, 0, NULL, NULL);
    //check_error(err, "clEnqueueReadBuffer C_row_ptr");

    err = clEnqueueReadBuffer(queue, C_col_idx_buf, CL_TRUE, 0, sizeof(C_col_idx), C_col_idx, 0, NULL, NULL);
    check_error(err, "clEnqueueReadBuffer C_col_idx");

    err = clEnqueueReadBuffer(queue, C_values_buf, CL_TRUE, 0, sizeof(C_values), C_values, 0, NULL, NULL);
    check_error(err, "clEnqueueReadBuffer C_values");

    //Prefix-sum to fit the C_row_ptr in the proper way
    //for (int i = 1; i <= A_nrows; i++)
      //  C_row_ptr[i] += C_row_ptr[i-1];
    
    // Print the result matrix in CSR format
    printf("Result Matrix C in CSR Format:\n");
    printf("C_row_ptr: ");
    for (int i = 0; i <= A_nrows; i++) {
        printf("%d ", C_row_ptr[i]);
    }
    printf("\n");

    printf("C_col_idx: ");
    for (int i = 0; i < C_row_ptr[A_nrows]; i++) {
        printf("%d ", C_col_idx[i]);
    }
    printf("\n");

    printf("C_values: ");
    for (int i = 0; i < C_row_ptr[A_nrows]; i++) {
        printf("%.2f ", C_values[i]);
    }
    printf("\n");

    // Clean up
    //clReleaseMemObject(A_row_ptr_buf);
    clReleaseMemObject(A_col_idx_buf);
    //clReleaseMemObject(A_values_buf);
    clReleaseMemObject(B_row_ptr_buf);
    clReleaseMemObject(B_col_idx_buf);
    clReleaseMemObject(B_values_buf);
    clReleaseMemObject(C_row_ptr_buf);
    clReleaseMemObject(C_col_idx_buf);
    clReleaseMemObject(C_values_buf);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);

    return 0;
}