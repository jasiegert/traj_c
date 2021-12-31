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
void linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R);

/**
* @brief Apply linear regression on data supplied as a 2D array with two columns.
*
* Applies a linear regression to the data points supplied as x and y.
* First column is assumed to be x, second column to be y.
* Only fit data in an interval defined by start_point and end_point.
*
* @param[in] n Number of points in the data.
* @param[in] ar array with data points in two columns (x and y).
* @param[in] start_point starting point expressed as fraction of total length, data points up until that point will be ignored.
* @param[in] end_point end point expressed as fraction of total length, data points beyond that point will be ignored.
* @param[out] slope slope of the linear fit.
* @param[out] intercept intercept with the y-axis of the linear fit.
* @param[out] R sum of all squared differences between data points and linear regression -> measure of 
*/
void linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept, float *R);

/**
* @brief Saves a float 2D array of arbitrary size into a csv-file.
*
* @param[in] outputname name of the file to be written (or replaced).
* @param[in] col_no number of columns.
* @param[in] row_no number of rows.
* @param[in] outputarray 2D array to be saved (or a pointer to its first element).
* @param[in] headerstring string to be printed into the first line of the csv or NULL if no header should be printed.
* @return 0 if the file could be written, 1 if it couldn't be opened.
*/
int savecsv(char *outputname, int col_no, int row_no, float outputarray[col_no][row_no], char *headerstring);

#endif /* MATHTOOLS_H */
