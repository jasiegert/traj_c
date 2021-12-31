/**
* @file oacf.h
* @author Johnny Alexander Jimenez Siegert
* @brief Functions to compute the orientational auto-correlation function (OACF) of a certain atom type in a trajectory.
*/

#ifndef OACF_H
#define OACF_H

/**
* @brief Calculates the orientational auto-correlation function (OACF) of the bonds between atoms of the specified types.
*
* Calculates the auto-correlation function of the bonds between the specified atom types.
* One bond will be drawn from each atom of type 1 to its closest neighbor of type 2.
*
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj atom coordinates of the trajectory.
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] pbc periodic boundary conditions as a 3x3 float array
* @param[in] target_atom_1 atomic number of the atom type, which will be the source of the bond.
* @param[in] target_atom_2 atomic number of the atom type, which will be the end of the bond.
* @param[in] resolution number of correlation times to sample.
* @param[in] timerange fraction of trajectory length, up to which correlation times will be generated.
* @param[out] oacf array holding correlation time and corresponding auto-correlation in each row.
* @param[out] output string to be printed and added to output file.
*/
int oacf_overall(
                int frame_no, \
                int atom_no, \
                float traj[frame_no][atom_no][3], \
                int atom[atom_no], \
                float pbc[3][3], \
                int target_atom_1, \
                int target_atom_2, \
                int resolution, \
                float timerange, \
                float oacf[resolution][2], \
                char *output\
                );

#endif /* OACF_H */
