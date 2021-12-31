/** @file mathtools.h
* @brief Function prototypes for simple mathematical functions.
*/


#ifndef MATHTOOLS_H
#define MATHTOOLS_H

/**
* @brief Apply linear regression on data supplied as two 1D float arrays.
*
* Applies a linear regression to the data points supplied as x and y.
* Only fit data in an interval defined by start_point and end_point.
*
* @param[in] n Number of points in the data.
* @param[in] x x-values of each data point.
* @param[in] y y-values of each data point.
* @param[in] start_point starting point expressed as fraction of total length, data points up until that point will be ignored.
* @param[in] end_point end point expressed as fraction of total length, data points beyond that point will be ignored.
* @param[out] slope slope of the linear fit.
* @param[out] intercept intercept with the y-axis of the linear fit.
* @param[out] R sum of all squared differences between data points and linear regression -> measure of 
*/
int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);
int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept, float *R);
int savecsv(char *outputname, int col_no, int row_no, float outputarray[col_no][row_no], char *headerstring);

#endif /* MATHTOOLS_H */
