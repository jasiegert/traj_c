/**
* @file rdf.h
* @author Johnny Alexander Jimenez Siegert
* @brief Functions to compute the radial distribution function (RDF) in a trajectory.
*/

#ifndef RDF_H
#define RDF_H

/**
* @brief Calculates the RDF from a trajectory.
*
* Calculates the radial-distribution function for the two atom types passed.
* Distances will be binned into the specified number of bins between min and max distance.
*
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj atom coordinates of the trajectory.
* @param[in] pbc periodic boundary conditions as a 3x3 float array
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] atom_1 first atom type of the RDF pair.
* @param[in] atom_2 second atom type of the RDF pair.
* @param[in] min_distance smallest distance to consider.
* @param[in] max_distance largest distance to consider.
* @param[in] bin number of bins and thus number of output values.
* @param[out] rdf sampled distances (average of each bin) and corresponding RDF.
* @param[out] output output string to be passed into stdout and output file.
*/
int rdf_overall(
                int frame_no, \ 
                int atom_no, \
                float traj[frame_no][atom_no][3], \
                float pbc[3][3], \
                int atom[atom_no], \
                int atom_1, \
                int atom_2, \
                float min_distance, \
                float max_distance, \
                int bin, \
                float rdf[bin][2], \
                char *output \
        );


/**
* @brief Calculates the intermolecular RDF from a trajectory.
*
* Calculates the radial-distribution function for the two atom types passed.
* Distances will be binned into the specified number of bins between min and max distance.
* Atoms will be assigned to molecules based on their closest atom of specified central atom type.
* Distances between atoms within the same molecule will be discarded.
*
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj atom coordinates of the trajectory.
* @param[in] pbc periodic boundary conditions as a 3x3 float array
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] atom_1 first atom type of the RDF pair.
* @param[in] atom_2 second atom type of the RDF pair.
* @param[in] atom_central central atom type used to assign atoms to molecules.
* @param[in] min_distance smallest distance to consider.
* @param[in] max_distance largest distance to consider.
* @param[in] bin number of bins and thus number of output values.
* @param[out] rdf sampled distances (average of each bin) and corresponding RDF.
* @param[out] output output string to be passed into stdout and output file.
*/
int rdf_intermolecular(
                int frame_no, \
                int atom_no, \
                float traj[frame_no][atom_no][3], \
                float pbc[3][3], \
                int atom[atom_no], \
                int atom_1, \
                int atom_2, \
                int atom_central, \
                float min_distance, \
                float max_distance, \
                int bin, \
                float rdf[bin][2], \
                char *output \
        );

#endif /* RDF_H */
