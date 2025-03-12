// Matrices A and B in CSR format
    /*
        [ 0  1  0  0 ]
        [ 1  0  0  0 ]
        [ 0  0  0  1 ]
        [ 0  0  1  0 ]
    */
    //Indexes (row):   0  1  2  3  4
    //int A_row_ptr[] = {0, 1, 2, 3, 4};
    int A_col_idx[] = {1, 0, 3, 2};
    //float A_values[] = {1.0, 1.0, 1.0, 1.0};
    int A_nrows = 4;
    int A_ncols = 4;

    /*      .  .
       -> [ 1  2  0  0 ]
        [ 2  1  0  0 ]
        [ 0  0  3  4 ]
        [ 0  0  4  3 ]
    */
                    // 0
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

//Compulação: A * B = C
[ 0  1  0  0 ]
[ 1  0  0  0 ]
[ 0  0  0  1 ] 
[ 0  0  1  0 ]
int A_col_idx[] = {1, 0, 3, 2};
*
[ 1  2  0  0 ]
[ 2  1  0  0 ]
[ 0  0  3  4 ]
[ 0  0  4  3 ]
int B_row_ptr[] = {0, 2, 4, 6, 8};
int B_col_idx[] = {0, 1, 0, 1, 2, 3, 2, 3};
float B_values[] = {1.0, 2.0, 2.0, 1.0, 3.0, 4.0, 4.0, 3.0};
=
[ 2  1  0  0 ]
[ 1  2  0  0 ]
[ 0  0  4  3 ]
[ 0  0  3  4 ]
int C_row_ptr[A_nrows + 1];
int C_col_idx[1024];
float C_values[1024];

//Computing thread: 0   1
int A_col_idx[] = {'1', *0*, 3, 2};

//Read
int B_row_ptr[] = {*0*, '2', 4, 6, 8};
int B_col_idx[] = {*0, 1*, '0, 1', 2, 3, 2, 3};
float B_values[] = {*1.0, 2.0*, '2.0, 1.0', 3.0, 4.0, 4.0, 3.0};

int C_row_ptr[] = {'0', *2*, 4, 6, 8}; //precomputed using prefix sum algorithm

//Write
int C_col_idx[] = {'0, 1', *0, 1*, 2, 3, 2, 3};
float C_values[] = {'2.0, 1.0', *1.0, 2.0*, ?, ?, ?, ?};



- Multiplicação de matrixes
- C = P X D
> Melhor O(nnz(D)) -> extremamente esparsa
> Pior O(P.rows x D.cols) -> Densa
>> multiplicação padrão é n^3