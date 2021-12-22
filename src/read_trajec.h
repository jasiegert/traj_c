#ifndef READTRAJEC_H
#define READTRAJEC_H

/**
* Counts number of atoms in an xyz-file.
* @param[in] name path to the xyz-file.
* @return number of atoms in the trajectory.
*/
int get_atoms(char *name);

/**
* Counts number of lines in a text file.
* @param[in] name path to the file.
* @return number of lines in the file.
*/
int get_lines(char *name);
int readxyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);
int readpbc(char *name, float pbc[3][3]);
int removecom(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float trajcom[frame_no][atom_no][3]);
int printarray(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);
int writexyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);
int savecsv(char *outputname, int col_no, int row_no, float outputarray[col_no][row_no]);

#endif /* READTRAJEC_H */
