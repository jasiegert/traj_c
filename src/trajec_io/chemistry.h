/**
* @file chemistry.h
* @brief Functions to get information on elements and do simple trajectory calculations.
* @ingroup TrajectoryIO
* @author Johnny Alexander Jimenez Siegert
*/

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

typedef struct atominfo
{
    int atom_no;
    float atom_mass;
    char symbol[3];
} atominfo;

/**
* @brief Fills a preallocated atominfo struct with information corresponding to a supplied atomic number.
* @param[in] atom_no atomic number
* @param[out] this_atom atominfo struct containing atomic number, atomic mass and element symbol
* @return 0 if atom number was found in range 1 to 118, otherwise 1
*/
int no_to_atominfo(int atom_no, atominfo *this_atom);

/**
* @brief Returns the atomic number for a given element name.
* @param[in] element string containing the element symbol (capitalized)
* @return (int) atomic number
*/
int element_to_no(char *element);

/**
* @brief Finds the element symbol corresponding to an atomic number.
* @param[in] atom_no atomic number.
* @param[out] element string containing the element symbol (capitalized)
*/
int no_to_element(int atom_no, char* element);

/**
* @brief Returns the atomic mass for a given element.
* @param[in] atom_no atomic number of the element.
* @return float atomic mass of the element (mass of an atom in atomic mass unit u).
*/
float no_to_mass(int atom_no);

/**
* @brief Calculates the distance between two points in real space following the minimum image convention.
* @param[in] coord_1, coord_2 Cartesian coordinates of each point.
* @param[in] pbc periodic boundary conditions represented as 3 cell vectors stacked into a 3x3 matrix.
* @return (float) distance between the two points.
* @warning Currently only supports orthogonal pbc (ignores non-diagonal values).
* @todo Add support for non-orthogonal pbc (as in rewrite it using matrix inversion and multiplication).
*/
float pbc_dist(float coord_1[3], float coord_2[3], float pbc[3][3]);

/**
* @brief Calculates the distance between two points in real space following the minimum image convention in an orthogonal box.
* @param[in] coord_1, coord_2 Cartesian coordinates of each point.
* @param[in] pbc periodic boundary conditions represented as 3 cell vectors stacked into a 3x3 matrix, off-diagonal elements must be zero.
* @return (float) distance between the two points.
*/
float pbc_dist_orthogonal(float coord_1[3], float coord_2[3], float pbc[3][3]);

/**
* @brief Calculates the distance between two points in real space following the minimum image convention.
* @param[in] coord_1, coord_2 Cartesian coordinates of each point.
* @param[in] pbc periodic boundary conditions represented as 3 cell vectors stacked into a 3x3 matrix, off-diagonal elements must be zero.
* @return (float) distance between the two points.
*/
float pbc_dist_triclinic(float coord_1[3], float coord_2[3], float pbc[3][3]);

/**
* @brief Calculates the distance between two points in real space following the minimum image convention.
* @param[in] pbc periodic boundary conditions represented as 3 cell vectors stacked into a 3x3 matrix, off-diagonal elements must be zero.
* @return (float) maximum distance for pbc_dist_triclinic to be guaranteed to return correct results
*/
float max_distance_pbc_dist_triclinic(float pbc[3][3]);

/**
* @brief Finds the closest atom ("neighbor") to a given point in a frame.
*
* Given a point and a frame (list of atom coordinates), returns the index of the entry in frame with the lowest distance to the point.\n
* If atom_neighbor is not 0, then only neighbors with it as atomic number may be returned. If 0, then all atoms are eligible.\n
* Distance is determined using \ref pbc_dist (see limitations there).
*
* @param[in] coord coordinates of the point, whose nearest neighbor should be found.
* @param[in] atom_no number of atoms in the frame.
* @param[in] traj_frame frame from the trajectory.
* @param[in] atom atomic numbers of all atoms in frame.
* @param[in] atom_neighbor atomic number of neighbor, which should be returned. Set to 0 to accept all atom types.
* @param[in] pbc periodic boundary conditions.
* @return (int) index of the atom closest to the given point in traj_frame. -1 if none is found.
*/
int nextneighbor_in_traj(float coord[3], int atom_no, float traj_frame[atom_no][3], int atom[atom_no], int atom_neighbor, float pbc[3][3]);

#endif /* CHEMISTRY_H */
