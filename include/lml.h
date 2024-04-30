/*
 * Lightweight Matrix Library (LML) - Header File
 *
 * File:     lml.h
 * Author:   James Bray
 * Repo:     https://github.com/James-Bray19/Lightweight-Matrix-Library
 *
 * Implementation file for the LML (Matrix Library) containing
 * definitions of functions declared in the lml.h header file.
 */

#ifndef LIGHTWEIGHTMATRIXLIBRARY_H
#define LIGHTWEIGHTMATRIXLIBRARY_H

// --------------- Matrix Structure ---------------

// the Matrix structure used as the fundamental data type in the LML library
typedef struct {
    int rows; // number of rows in the matrix
    int cols; // number of columns in the matrix.
    double **data; // pointer to a 2D array of doubles, which stores the matrix data
} Matrix;

// --------------- Generating Matrices ---------------

/// @brief returns a matrix filled with zeros of size rows x cols
/// @param rows number of rows
/// @param cols number of columns
/// @return matrix filled with zeros
Matrix *zeros(int rows, int cols);

/// @brief returns a matrix filled with a specific value of size rows x cols 
/// @param rows number of rows
/// @param cols number of columns
/// @return matrix filled with specified value
Matrix *constants(int rows, int cols, int value);

/// @brief returns an identity matrix of size size x size
/// @param size size of the identity matrix
/// @return identity matrix
Matrix *identity(int size);

/// @brief returns a matrix filled with random values of size rows x cols
/// @param rows number of rows
/// @param cols number of columns
/// @return matrix filled with random values
Matrix *random(int rows, int cols);

/// @brief returns a matrix filled with random values of size rows x cols
/// @param rows number of rows
/// @param cols number of columns
/// @param lower lower bound for random values
/// @param upper upper bound for random values
/// @return matrix filled with random values
Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);

// --------------- Retrieving Data ---------------

/// @brief copy the input matrix
/// @param mat input matrix
/// @return copied matrix
Matrix *copy(Matrix *mat);

/// @brief returns the specified row of the matrix
/// @param mat input matrix
/// @param row row index
/// @return row of the matrix
Matrix *get_row(Matrix *mat, int row);

/// @brief returns the specified column of the matrix
/// @param mat input matrix
/// @param col column index
/// @return column of the matrix
Matrix *get_col(Matrix *mat, int col);

/// @brief returns the values in the lower triangular matrix
/// @param mat input matrix
/// @return lower triangular matrix
Matrix *get_lower(Matrix *mat);

/// @brief returns the values in the upper triangular matrix
/// @param mat input matrix
/// @return upper triangular matrix
Matrix *get_upper(Matrix *mat);

/// @brief returns the submatrix of the input matrix
/// @param mat input matrix
/// @param row starting row index
/// @param col starting column index
/// @param rows number of rows
/// @param cols number of columns
/// @return submatrix
Matrix *get_submatrix(Matrix *mat, int row, int col, int rows, int cols);

// --------------- Matrix Operations ---------------

/// @brief returns the determinant of the input matrix
/// @param mat input matrix
/// @return determinant of the matrix
double det(Matrix *mat);

/// @brief returns the transpose of the input matrix
/// @param mat input matrix
/// @return transposed matrix
Matrix *transpose(Matrix *mat);

/// @brief returns the result of matrix addition of two matrices
/// @param mat1 first matrix
/// @param mat2 second matrix
/// @return result of matrix addition
Matrix *add(Matrix *mat1, Matrix *mat2);

/// @brief returns the result of matrix multiplication of two matrices
/// @param mat1 first matrix
/// @param mat2 second matrix
/// @return result of matrix multiplication
Matrix *multiply(Matrix *mat1, Matrix *mat2);

/// @brief scales the matrix by multiplying each element by a scalar value
/// @param mat input matrix
/// @param scalar scalar value
/// @return scaled matrix
Matrix *scalar_multiply(Matrix *mat, double scalar);

/// @brief shifts the matrix by adding a scalar value to each element
/// @param mat input matrix
/// @param scalar scalar value
/// @return shifted matrix
Matrix *scalar_add(Matrix *mat, double scalar);

/// @brief returns the result of matrix multiplication of two matrices
/// @param mat input matrix
/// @param L pointer to store the lower triangular matrix
/// @param U pointer to store the upper triangular matrix
/// @return L and U matrices
void LU_decompose(Matrix *mat, Matrix **L, Matrix **U);

/// @brief returns the result of matrix multiplication of two matrices
/// @param mat1 first matrix
/// @param mat2 second matrix
/// @return result of matrix multiplication
Matrix *solve(Matrix *mat1, Matrix *mat2);

/// @brief returns the inverse of the input matrix
/// @param mat input matrix
/// @return inverse of the matrix
Matrix *inverse(Matrix *mat);

// --------------- In-Place Operations ---------------

/// @brief applies a function to each element of the matrix
/// @param mat input matrix
/// @param function function to apply
/// @return matrix with function applied
void map(Matrix *mat, double (*function)(double));

/// @brief sets the values of a specific row in the matrix
/// @param mat input matrix
/// @param row_index row index
/// @param row_values matrix of row values
void set_row(Matrix *mat, int row_index, Matrix *row_values);

/// @brief sets the values of a specific column in the matrix
/// @param mat input matrix
/// @param col_index column index
/// @param col_values matrix of column values
void set_col(Matrix *mat, int col_index, Matrix *col_values);

/// @brief sets the values of a submatrix inside a larger matrix
/// @param mat input matrix
/// @param row starting row index
/// @param col starting column index
/// @param sub submatrix to insert
void set_submatrix(Matrix *mat, int row, int col, Matrix *sub);

/// @brief removes a row from the matrix
/// @param mat input matrix
/// @param row row index
void remove_row(Matrix *mat, int row);

/// @brief removes a column from the matrix
/// @param mat input matrix
/// @param col column index
void remove_col(Matrix *mat, int col);

/// @brief inserts a new row with the provided values at the specified row index
/// @param mat destination matrix
/// @param row_index row index
/// @param row_values matrix of row values
void insert_row(Matrix *mat, int row, Matrix *row_values);

/// @brief inserts a new column with the provided values at the specified column index
/// @param mat destination matrix
/// @param col_index column index
/// @param col_values matrix of column values
void insert_col(Matrix *mat, int col, Matrix *col_values);

/// @brief appends rows from mat2 to mat1
/// @param mat1 destination matrix
/// @param mat2 source matrix
void append_rows(Matrix *mat1, Matrix *mat2);

/// @brief appends columns from mat2 to mat1
/// @param mat1 destination matrix
/// @param mat2 source matrix
void append_cols(Matrix *mat1, Matrix *mat2);

// --------------- Miscellaneous Functions ---------------

/// @brief displays the matrix
/// @param mat input matrix
void display(Matrix *mat);

/// @brief releases the memory allocated for the matrix
/// @param mat input matrix
void release(Matrix *mat);

#endif /* LIGHTWEIGHTMATRIXLIBRARY_H */