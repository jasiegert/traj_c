/** @file mathtools.h
* @brief Function prototypes for simple mathematical functions.
*/


#ifndef MATHTOOLS_H
#define MATHTOOLS_H

int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);
int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept, float *R);
int savecsv(char *outputname, int col_no, int row_no, float outputarray[col_no][row_no]);

#endif /* MATHTOOLS_H */
