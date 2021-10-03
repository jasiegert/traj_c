#ifndef CHEMISTRY_H
#define CHEMISTRY_H

int element_to_no(char *element);
int no_to_element(int atom_no, char* element);
float no_to_mass(int atom_no);
int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);
int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept,     float *R);

#endif /* CHEMISTRY_H */
