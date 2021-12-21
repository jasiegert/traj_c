#ifndef OACF_H
#define OACF_H

int oacf_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float pbc[3][3], int target_atom_1, int target_atom_2, int resolution, float timerange, float oacf[resolution][2], char *output);

#endif /* OACF_H */
