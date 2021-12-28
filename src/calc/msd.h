/**
* @file msd.h
* @author Johnny Alexander Jimenez Siegert
* @brief Functions to compute the mean square displacement (MSD) of a certain atom type in a trajectory.
*/

#ifndef MSD_H
#define MSD_H

/**
* @brief Calculates the mean square displacement (MSD) of a certain atom type in a trajectory.
*
* 
*
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj atom coordinates of the trajectory.
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] target_atom atomic number of the atom type, whose MSD should be calculated.
* @param[in] resolution number of correlation times to sample.
* @param[in] timerange fraction of trajectory length, up to which correlation times will be generated.
* @param[out] msd array holding correlation time and corresponding MSD in each row.
* @param[out] output string to be printed and added to output file.
*/
int msd_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, int resolution, float timerange, float msd[resolution][2], char *output);
/**
* @brief Calculates the mean square displacement (MSD) of a certain atom type in a trajectory.
*
* 
*
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj atom coordinates of the trajectory.
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] target_atom atomic number of the atom type, whose MSD should be calculated.
* @param[in] timerange fraction of trajectory length, up to which correlation times will be generated.
* @param[out] msd array holding correlation time and corresponding MSD in each row.
* @param[out] output string to be printed and added to output file.
*/
int msd_fft(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, float timerange, float msd[frame_no][2], char *output);

#endif /* MSD_H */
