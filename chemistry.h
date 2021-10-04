#ifndef CHEMISTRY_H
#define CHEMISTRY_H

int element_to_no(char *element);
int no_to_element(int atom_no, char* element);
float no_to_mass(int atom_no);
int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);
int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept,     float *R);
float pbc_dist(float coord_1[3], float coord_2[3], float pbc[3][3]);
int pbc_coord(float coord[3], float coord_pbc[3], float pbc[3][3]);
int pbc_traj(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float traj_pbc[frame_no][atom_no][3], float pbc[3][3]);

#endif /* CHEMISTRY_H */
