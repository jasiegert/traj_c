/** @file read_trajec.h
* @author Johnny Alexander Jimenez Siegert
* @brief Functions to read/write trajectories from/to xyz-files.
* @ingroup TrajectoryIO 
* @warning Currently only supports xyz-files. Other common file formats might be added later on.
*/

#ifndef READTRAJEC_H
#define READTRAJEC_H

/**
* @brief Counts number of atoms and frames in an xyz-file.
*
* Each frame should consist of 2 comment lines followed by atom_no entries.\n
* Number of atoms will be counted in the first frame. Frame number will be determined by counting the lines and dividing by line number per frame (atom_no + 2).
*
* @param[in] name path to the xyz-file.
* @param[out] atom_no number of atoms.
* @param[out] frame_no number of frames.
* @return 1 if something went wrong (e.g. file can't be opened or is not properly formatted), otherwise 0.
*/
int get_atom_and_frame_no(char *name, int *atom_no, int *frame_no);

void skipline(FILE *f);

/**
* @brief Reads content of an xyz-file or dat-file into a 3D-array.
*
* Reads a trajectory from an xyz-file. If an auxiliary dat-file (with .xyz replaced by .dat) exists, it will be used instead, otherwise it will be generated.\n
* Since atom is a 1D-array the pointer *atom_pointer can be used as is. For the trajectory casting to a 3D-array might be useful:\n
* \code
* float *trajectory_pointer;
* int frame_no, atom_no;
* int *atom;
* readtraj(char *name, &frame_no, &atom_no, &trajectory_pointer, &atom);
* float (*traj)[atom_no][3] = ( float (*)[atom_no][3] ) trajectory_pointer;
* \endcode
*
* @warning Atom (*atom_pointer) and trajectory (*trajectory_pointer) will be malloc'd, remember to free them!
* @param[in] name path to the xyz-file.
* @param[in] frame_no_pointer number of frames in the trajectory.
* @param[in] atom_no_pointer number of atoms in the trajectory
* @param[out] trajectory_pointer 3D-array to hold the trajectory coordinates.
* @param[out] atom_pointer atomic numbers of atoms in the trajectory in the same order as they appear in the first frame.
* @return int 1 if a problem occured, otherwise 0.
*/
int readtraj(char *name, int *frame_no_pointer, int *atom_no_pointer, float **trajectory_pointer, int **atom_pointer);

/**
* @brief Reads content of an xyz-file into a 3D-array.
*
* Reads xyz-file. For usage see \ref readtraj.
*
* @param[in] name path to the xyz-file.
* @param[in] frame_no_pointer number of frames in the trajectory.
* @param[in] atom_no_pointer number of atoms in the trajectory
* @param[out] trajectory_pointer 3D-array to hold the trajectory coordinates.
* @param[out] atom_pointer atomic numbers of atoms in the trajectory in the same order as they appear in the first frame.
* @return int 1 if a problem occured, otherwise 0.
*/
int readxyz(char *name, int *frame_no_pointer, int *atom_no_pointer, float **trajectory_pointer, int **atom_pointer);

/**
* @brief Reads content of an auxiliary dat-file into a 3D-array.
*
* Reads dat-file. For usage see \ref readtraj.
*
* @param[in] datname path to the xyz-file.
* @param[in] frame_no_pointer number of frames in the trajectory.
* @param[in] atom_no_pointer number of atoms in the trajectory
* @param[out] trajectory_pointer 3D-array to hold the trajectory coordinates.
* @param[out] atom_pointer atomic numbers of atoms in the trajectory in the same order as they appear in the first frame.
* @return int 1 if a problem occured, otherwise 0.
*/
int readdat(char *datname, int *frame_no_pointer, int *atom_no_pointer, float **trajectory_pointer, int **atom_pointer);

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
