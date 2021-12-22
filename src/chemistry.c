#include <string.h>
#include <stdio.h>
#include <math.h>

int element_to_no(char *element)
{
    int atom_no;
    if (strcmp(element, "H") == 0) atom_no = 1;
    else if (strcmp(element, "Li") == 0) atom_no = 3;
    else if (strcmp(element, "O") == 0) atom_no = 8;
    else if (strcmp(element, "C") == 0) atom_no = 6;
    else if (strcmp(element, "N") == 0) atom_no = 7;
    else if (strcmp(element, "O") == 0) atom_no = 8;
    else if (strcmp(element, "P") == 0) atom_no = 15;
    else if (strcmp(element, "S") == 0) atom_no = 16;
    else if (strcmp(element, "Cl") == 0) atom_no = 17;
    else if (strcmp(element, "Cs") == 0) atom_no = 55;
    else return -1;

    return atom_no;
}

int no_to_element(int atom_no, char* element)
{
    switch (atom_no)
    {
        case 1: strcpy(element, "H"); break;
        case 3: strcpy(element, "Li"); break;
        case 6: strcpy(element, "C"); break;
        case 7: strcpy(element, "N"); break;
        case 8: strcpy(element, "O"); break;
        case 15: strcpy(element, "P"); break;
        case 16: strcpy(element, "S"); break;
        case 17: strcpy(element, "Cl"); break;
        case 55: strcpy(element, "Cs"); break;
        default: return -1;
    }
    return 0;
}

float no_to_mass(int atom_no)
{
    float atom_mass;
    switch (atom_no)
    {
        case 1: atom_mass = 1.0; break;
        case 3: atom_mass = 0.0; break;
        case 6: atom_mass = 12.0; break;
        case 7: atom_mass = 14.0; break;
        case 8: atom_mass = 16.0; break;
        case 15: atom_mass = 31.0; break;
        case 16: atom_mass = 32.0; break;
        case 17: atom_mass = 35.45; break;
        case 55: atom_mass = 133.0; break;
        default: return -1.0;
    }
//    if (atom_mass == 0)
//    {
//        char atom_label[3];
//        no_to_element(atom_no, atom_label);
//        printf("Warning: atom mass of %s is %f.", atom_label, atom_mass);
 //   }
    return atom_mass;
}

int linregress(int n, float x[n], float y[n], float start_point, float end_point, float *slope, float *intercept, float *R)
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

    return 0;
}

int linregress_array(int n, float ar[n][2], float start_point, float end_point, float *slope, float *intercept, float *R)
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

    return 0;
}

float pbc_dist(float coord_1[3], float coord_2[3], float pbc[3][3])
{
    float diff[3];
    float dist = 0;
    for (int i = 0; i < 3; i++)
    {
        diff[i] = coord_2[i] - coord_1[i];
        if (diff[i] < 0) 
            diff[i] *= -1;
        diff[i] = diff[i] - floor(diff[i] / pbc[i][i]) * pbc[i][i];
        if (diff[i] > pbc[i][i]/2)
        {
            diff[i] = pbc[i][i] - diff[i];
        }

        dist += diff[i] * diff[i];
    }
    dist = sqrt(dist);
    return dist;
}

// Finds closest neighbor of type atom_neighbor in trajectory defined by traj_frame and atom to coordinates coord; returns -1 if nothing is found
int nextneighbor_in_traj(float coord[3], int atom_no, float traj_frame[atom_no][3], int atom[atom_no], int atom_neighbor, float pbc[3][3])
{
    int index_neighbor = -1;
    float min_distance = 10.0; //sqrt(sqr(pbc[0][0]) + sqr(pbc[1][1]) + sqr(pbc[2][2])); // longest distance possible in orthogonal cell with side lengths pbc[i][i]
    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_neighbor)
        {
            float dist = pbc_dist(coord, traj_frame[i], pbc);
            if (dist < min_distance)
            {
                min_distance = dist;
                index_neighbor = i;
            }
        }
        else 
            continue;
    }

    return index_neighbor;
}

int pbc_coord(float coord[3], float coord_pbc[3], float pbc[3][3])
{
    return 0;
}

int pbc_traj(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float traj_pbc[frame_no][atom_no][3], float pbc[3][3])
{
    return 0;
}
