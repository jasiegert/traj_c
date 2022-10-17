/**
* @file matrix_and_vectors.h
* @brief Selected operations on 3x3 matrices and 3 vectors, which are necessary to transform real to fractional coordinates.
* @ingroup TrajectoryIO
* @author Johnny Alexander Jimenez Siegert
*/

#ifndef MATRICES_H
#define MATRICES_H

/**
* @brief Calculates cofactor matrix.
* @param[in] mat original matrix (3x3)
* @param[out] cof cofactor matrix (3x3)
*/
void matrix33_cofactors(float mat[3][3], float cof[3][3]);
/**
* @brief Creates transposed matrix.
* @param[in] mat original matrix (3x3)
* @param[out] transpose transposed matrix (3x3)
*/
void matrix33_transpose(float mat[3][3], float transpose[3][3]);
/**
* @brief Calculcate product of two matrices c = a * b.
* @param[in] a first matrix (3x3)
* @param[in] b second matrix (3x3)
* @param[out] c product matrix (3x3)
*/
void matrix33_multiplication(float a[3][3], float b[3][3], float c[3][3]);
/**
* @brief Calculcates inverse of a matrix.
* @param[in] mat original matrix (3x3)
* @param[out] inv inverse matrix (3x3)
*/
void matrix33_inverse(float mat[3][3], float inv[3][3]);
/**
* @brief Calculcate product vector of matrix and vector.
* @param[in] mat matrix (3x3)
* @param[in] vec vector (3)
* @param[out] out product vector mat * vec
*/
void matrix33_vector3_multiplication(float mat[3][3], float vec[3], float out[3]);

#endif /* MATRICES_H */
