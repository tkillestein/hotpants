/**
 * @file allocate.c
 * @brief Implementation of safe memory allocation wrappers.
 *
 * Provides xmalloc/xcalloc for consistent error handling and
 * contiguous multi-dimensional array allocation for efficiency.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allocate.h"

/**
 * Allocate memory with automatic error checking and diagnostic message.
 * Exits the program if allocation fails.
 */
void *xmalloc(size_t size) {
    void *ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "[ERROR] Memory allocation failed: malloc(%zu bytes) returned NULL\n", size);
        exit(1);
    }
    return ptr;
}

/**
 * Allocate and zero-initialize memory with automatic error checking.
 * Exits the program if allocation fails.
 */
void *xcalloc(size_t count, size_t size) {
    void *ptr = calloc(count, size);
    if (!ptr) {
        fprintf(stderr, "[ERROR] Memory allocation failed: calloc(%zu, %zu) returned NULL\n", count, size);
        exit(1);
    }
    return ptr;
}

/**
 * Allocate a 2D matrix as a single contiguous block of memory.
 *
 * The matrix is stored row-major in a single malloc, which improves cache
 * locality and simplifies cleanup. The returned pointer is an array of
 * row pointers that index into the contiguous data block.
 *
 * Layout:
 *   matrix[0] -> points to row 0 (data[0*cols..0*cols+cols-1])
 *   matrix[1] -> points to row 1 (data[1*cols..1*cols+cols-1])
 *   etc.
 *
 * Free with: free_matrix_contiguous(matrix)
 */
double **alloc_matrix_contiguous(int rows, int cols) {
    double *data;
    double **ptrs;
    int i;

    /* Allocate the flat data array */
    data = (double *)xcalloc((size_t)rows * cols, sizeof(double));

    /* Allocate the array of row pointers */
    ptrs = (double **)xmalloc((size_t)rows * sizeof(double *));

    /* Initialize each row pointer to point into the data block */
    for (i = 0; i < rows; i++) {
        ptrs[i] = data + (size_t)i * cols;
    }

    return ptrs;
}

/**
 * Free a 2D matrix allocated by alloc_matrix_contiguous().
 *
 * Only two free() calls are needed: one for the row pointer array,
 * one for the contiguous data block (stored as ptrs[0]).
 */
void free_matrix_contiguous(double **matrix) {
    if (!matrix) return;

    /* Free the contiguous data block (pointed to by matrix[0]) */
    free(matrix[0]);

    /* Free the row pointer array */
    free(matrix);
}

/**
 * Allocate a contiguous array of vectors (e.g., for stamp->vectors).
 *
 * Each vector is stored contiguously in a single data block, accessed via
 * an array of pointers. This is more efficient than nVectors separate
 * malloc() calls.
 *
 * Example usage (stamp allocation):
 *   stamp->vectors = alloc_vector_array(nCompKer + nBGVectors, fwKSStamp*fwKSStamp);
 *
 * Layout:
 *   vectors[0] -> data[0*vectorSize .. 0*vectorSize+vectorSize-1]
 *   vectors[1] -> data[1*vectorSize .. 1*vectorSize+vectorSize-1]
 *   etc.
 *
 * Free with: free_vector_array(vectors)
 */
double **alloc_vector_array(int nVectors, int vectorSize) {
    double *data;
    double **ptrs;
    int i;

    /* Allocate the contiguous data block for all vectors */
    data = (double *)xcalloc((size_t)nVectors * vectorSize, sizeof(double));

    /* Allocate the array of vector pointers */
    ptrs = (double **)xmalloc((size_t)nVectors * sizeof(double *));

    /* Initialize each vector pointer to point into the data block */
    for (i = 0; i < nVectors; i++) {
        ptrs[i] = data + (size_t)i * vectorSize;
    }

    return ptrs;
}

/**
 * Free a vector array allocated by alloc_vector_array().
 *
 * Only two free() calls: one for the pointer array, one for the
 * contiguous data block.
 */
void free_vector_array(double **vectors) {
    if (!vectors) return;

    /* Free the contiguous data block */
    free(vectors[0]);

    /* Free the pointer array */
    free(vectors);
}
