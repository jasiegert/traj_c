/**
* @file matrix_and_vectors.h
* @brief Selected operations on 3x3 matrices and 3 vectors.
* @ingroup TrajectoryIO
* @author Johnny Alexander Jimenez Siegert
*/

#ifndef MATRICES_H
#define MATRICES_H

/**
* @brief Finds cofactor matrix to a given 3x3 matrix.
* @param[in] mat original matrix
* @param[out] cofactors cofactor matrix
*/
void matrix33_adjoint(float mat[3][3], float adj[3][3]);
void matrix33_transpose(float mat[3][3], float transpose[3][3]);
void matrix33_multiplication(float a[3][3], float b[3][3], float c[3][3]);
void matrix33_inverse(float mat[3][3], float inv[3][3]);
void matrix33_vector3_multiplication(float mat[3][3], float vec[3], float out[3]);

#endif /* MATRICES_H */
