#ifndef READTRAJEC_H
#define READTRAJEC_H

int get_atoms(char *name);
int get_lines(char *name);
int readxyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);
int readpbc(char *name, float pbc[3][3]);
int removecom(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float trajcom[frame_no][atom_no][3]);
int printarray(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);
int writexyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no]);

#endif /* READTRAJEC_H */
