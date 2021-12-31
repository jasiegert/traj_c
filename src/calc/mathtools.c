#include <stdio.h>
#include <math.h>

void linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R)
{
    float xsum = 0, x2sum = 0, ysum = 0, y2sum = 0, xysum = 0;
    int start_i = round(start_point * n);
    int end_i = round(end_point * n);
    int diff_i = end_i - start_i;
    for (int i = start_i; i < end_i; i++)
    {
        xsum += x[i];
        x2sum += x[i] * x[i];
        ysum += y[i];
        y2sum += y[i] * y[i];
        xysum += x[i] * y[i];
    }

    *slope = (diff_i * xysum - xsum * ysum) / (diff_i * x2sum - xsum * xsum);
    *intercept = (x2sum * ysum - xsum * xysum) / (diff_i * x2sum - xsum * xsum);
    *R = (diff_i * xysum - xsum * ysum) / (sqrt(diff_i * x2sum - xsum * xsum) * sqrt(diff_i * y2sum - ysum * ysum));
}

void linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept, float *R)
{
    float xsum = 0.0, x2sum = 0.0, ysum = 0.0, y2sum = 0.0, xysum = 0.0;
    int start_i = round(start_point * n);
    int end_i = round(end_point * n);
    int diff_i = end_i - start_i;
    for (int i = start_i; i < end_i; i++)
    {
        xsum += ar[i][0];
        x2sum += ar[i][0] * ar[i][0];
        ysum += ar[i][1];
        y2sum += ar[i][1] * ar[i][1];
        xysum += ar[i][0] * ar[i][1];
    }

    *slope = (diff_i * xysum - xsum * ysum) / (diff_i * x2sum - xsum * xsum);
    *intercept = (x2sum * ysum - xsum * xysum) / (diff_i * x2sum - xsum * xsum);
    *R = (diff_i * xysum - xsum * ysum) / (sqrt(diff_i * x2sum - xsum * xsum) * sqrt(diff_i * y2sum - ysum * ysum));
}

int savecsv(char *outputname, int col_no, int row_no, float outputarray[col_no][row_no], char *headerstring)
{
    FILE *output = fopen(outputname, "w");
    if (output == NULL)
    {
        return 1;
    }

    if (headerstring != NULL)
    {
        fputs(headerstring, output);
        fputs("\n", output);
    }

    for (int i = 0; i < col_no; i++)
    {
        for (int j = 0; j < row_no; j++)
        {
            fprintf(output, "%.10e    ", outputarray[i][j]);
        }
        fprintf(output, "\n");
    }

    fclose(output);
    return 0;
}
