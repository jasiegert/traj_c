/** @file read_trajec.h
* @author Johnny Alexander Jimenez Siegert
* @brief Functions to read/write trajectories from/to xyz-files.
* @ingroup TrajectoryIO 
* @warning Currently only supports xyz-files. Other common file formats might be added later on.
*/

#ifndef READTRAJEC_H
#define READTRAJEC_H

/**
* @brief Counts number of atoms in an xyz-file.
*
* @param[in] name path to the xyz-file.
* @return number of atoms in the trajectory.
*/
int get_atoms(char *name);

/**
* @brief Counts number of lines in a text file.
*
* @param[in] name path to the file.
* @return number of lines in the file.
*/
int get_lines(char *name);

int get_atom_and_frame_no(char *name, int *atom_no, int *frame_no);

void skipline(FILE *f);

int readtraj(char *name, int *frame_no, int *atom_no, float **trajectory_pointer, int **atom_pointer);

/**
* @brief Reads content of an xyz-file into a 3D-array.
*
* @param[in] name path to the xyz-file.
* @param[in] frame_no number of frames in the trajectory.
* @param[in] atom_no number of atoms in the trajectory
* @param[out] traj 3D-array to hold the trajectory coordinates.
* @param[out] atom atomic numbers of atoms in the trajectory in the same order as they appear in the first frame.
* @return int 1 if a problem occured, otherwise 0.
*/
int readxyz(char *name, int *frame_no_pointer, int *atom_no_pointer, float** trajectory_pointer, int **atom_pointer);

int readdat(char *datname, int *frame_no, int *atom_no, float **trajectory_pointer, int **atom_pointer);

/**
* @brief Reads file containing periodic boundary conditions into array.
* 
* Periodic boundary conditions have to be represented as three vectors.
* The file passed to this function should contain one vector represented as three floats on each line.
* Additional lines will be ignored.
*
* @param[in] name path to the file.
* @param[out] pbc 3x3 array containing each pbc vector as a row.
* @return 1 if a problem occured, otherwise 0.
*/
int readpbc(char *name, float pbc[3][3]);

/**
* @brief Removes center of mass movement from a trajectory
*
* Removes center of mass movement from the trajectory by setting the center of mass to (0, 0, 0) in each frame.
* This is necessary to correctly calculate the mean-square displacement.\n
* Source (traj) and target (trajcom) can be the same.
*
* @param[in] frame_no number of frames in the trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj coordinates of all atoms.
* @param[in] atom atomic numbers of all atoms.
* @param[out] trajcom coordinates of all atoms with the center of mass set to (0, 0, 0) in each frame.
*/
int removecom(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float trajcom[frame_no][atom_no][3]);

/**
* @brief Writes trajectory (atom labels and coordinates) into an xyz-file.
*
* @param[in] name path to the file, which should be written.
* @param[in] frame_no number of frames in the trajectory.
* @param[in] atom_no number of atoms in the trajectory.
* @param[in] traj coordinates.
* @param[in] atom atomic numbers.
*/
int writexyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);

#endif /* READTRAJEC_H */
