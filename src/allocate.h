/**
 * @file allocate.h
 * @brief Safe memory allocation wrappers with automatic error handling and
 *        optimized multi-dimensional array allocation.
 *
 * This module provides convenient wrappers around malloc/calloc that:
 * 1. Check for allocation failure and exit cleanly with diagnostic message
 * 2. Allocate contiguous 2D/3D arrays for better cache locality
 * 3. Handle cleanup with a single function call
 *
 * Usage: Replace direct malloc() calls with xmalloc(), and nested 2D array
 *        allocations with alloc_matrix_contiguous() for improved efficiency
 *        and reduced fragmentation.
 *
 * Example:
 *   double *flat = xcalloc(100*100, sizeof(double));
 *   double **matrix = alloc_matrix_contiguous(100, 100);
 */

#ifndef ALLOCATE_H
#define ALLOCATE_H

#include <stddef.h>

/* =====================================================================
 * Basic Safe Allocation Wrappers
 * =====================================================================
 * These functions wrap malloc/calloc and exit with a diagnostic message
 * on allocation failure, eliminating the need for repeated error checks.
 */

/**
 * Allocate memory like malloc(), but exit with error message on failure.
 *
 * @param size Number of bytes to allocate
 * @return Pointer to allocated memory (never NULL; exits on failure)
 */
void *xmalloc(size_t size);

/**
 * Allocate memory like calloc(), but exit with error message on failure.
 * Memory is zero-initialized.
 *
 * @param count Number of elements
 * @param size  Size of each element in bytes
 * @return Pointer to allocated, zero-initialized memory (never NULL)
 */
void *xcalloc(size_t count, size_t size);

/* =====================================================================
 * Contiguous Multi-Dimensional Array Allocation
 * =====================================================================
 * Allocate 2D and 3D arrays as single contiguous blocks, avoiding the
 * fragmentation and poor cache behavior of nested malloc loops.
 *
 * Example:
 *   // Old way: 3 mallocs, poor cache locality
 *   double **matrix = malloc(rows * sizeof(double*));
 *   for(i=0; i<rows; i++) matrix[i] = malloc(cols*sizeof(double));
 *
 *   // New way: 2 mallocs, contiguous data
 *   double **matrix = alloc_matrix_contiguous(rows, cols);
 */

/**
 * Allocate a 2D array in contiguous memory.
 *
 * The array data is stored as a single contiguous block for better cache
 * locality. The row pointers (matrix[i]) point into different sections
 * of this block.
 *
 * Memory layout: All elements in a single malloc, indexed as:
 *   matrix[i][j] = data[i*cols + j]
 *
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Pointer to array of row pointers; free with free_matrix_contiguous()
 */
double **alloc_matrix_contiguous(int rows, int cols);

/**
 * Free a 2D array allocated by alloc_matrix_contiguous().
 *
 * @param matrix Pointer returned from alloc_matrix_contiguous()
 */
void free_matrix_contiguous(double **matrix);

/**
 * Allocate a contiguous array of vectors (array of pointers to arrays).
 *
 * Used for stamp->vectors: allocates nVectors separate vector arrays,
 * each of size vectorSize, in a single contiguous block. Returns an array
 * of pointers to the start of each vector.
 *
 * Memory layout:
 *   vectors[i] points to element i*vectorSize in the data block
 *
 * @param nVectors   Number of vectors
 * @param vectorSize Size of each vector (in elements, not bytes)
 * @return Array of nVectors pointers; free with free_vector_array()
 */
double **alloc_vector_array(int nVectors, int vectorSize);

/**
 * Free a vector array allocated by alloc_vector_array().
 *
 * @param vectors Pointer returned from alloc_vector_array()
 */
void free_vector_array(double **vectors);

#endif  /* ALLOCATE_H */
