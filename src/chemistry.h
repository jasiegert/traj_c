#ifndef CHEMISTRY_H
#define CHEMISTRY_H


/**
* Returns the atomic number for a given element name.
* @param[in] element string containing the element symbol (capitalized)
* @return (int) atomic number
*/
int element_to_no(char *element);

/**
* Finds the element symbol corresponding to an atomic number.
* @param[in] atom_no atomic number.
* @param[out] element string containing the element symbol (capitalized)
*/
int no_to_element(int atom_no, char* element);


float no_to_mass(int atom_no);
int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);
int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept,     float *R);

/**
* Calculates the distance between two points in real space following the minimum image convention.
* @param[in] coord_1, coord_2 Cartesian coordinates of each point.
* @param[in] pbc periodic boundary conditions represented as 3 cell vectors stacked into a 3x3 matrix.
* @return (float) distance between the two points.
* @warning Currently only supports orthogonal pbc (ignores non-diagonal values).
* @todo Add support for non-orthogonal pbc (as in rewrite it using matrix inversion and multiplication).
*/
float pbc_dist(float coord_1[3], float coord_2[3], float pbc[3][3]);

int nextneighbor_in_traj(float coord[3], int atom_no, float traj_frame[atom_no][3], int atom[atom_no], int atom_neighbor, float pbc[3][3]);
int pbc_coord(float coord[3], float coord_pbc[3], float pbc[3][3]);
int pbc_traj(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float traj_pbc[frame_no][atom_no][3], float pbc[3][3]);

#endif /* CHEMISTRY_H */
