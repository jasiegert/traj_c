#include <stdio.h>
#include <math.h>

#include "../chemistry.h"

// Calculates orientational autocorrelation function of vector connecting each target_atom_1 and its closest target_atom_2 neighbor
int oacf_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float pbc[3][3], int target_atom_1, int target_atom_2, int resolution, float timerange, float oacf[resolution][2], char *output)
{
    float timestep = 0.5; // in fs

    // Count number of target atoms in trajectory and save indices of target atom and its closest target_atom_2
    int atom_1_no = 0;
    int atom_1_index[atom_no];
    int atom_2_index[atom_no];
    for (int j = 0; j < atom_no; j++)
    {
        if (atom[j] == target_atom_1)
        {
            atom_1_index[atom_1_no] = j;
            atom_2_index[atom_1_no] = nextneighbor_in_traj(traj[0][j], atom_no, traj[0], atom, target_atom_2, pbc);
            atom_1_no++;
        }
    }

    // Create correlation time for each resolution point 
    for (int i = 0; i < resolution; i++)
    {
        int corr_time = round((i + 1) / (float)resolution * timerange * frame_no);
        oacf[i][0] = corr_time * timestep / 1000;   // Correlation time in ps
        oacf[i][1] = 0;             // MSD for correlation time in pm^2

        // Iterate through each target_atom_1 and calculate oacf
        for (int j = 0; j < atom_1_no; j++)
        {
            int atom_1_tmp = atom_1_index[j];
            int atom_2_tmp = atom_2_index[j];

            // Iterate through correlation windows (shifted by one timestep each)
            for (int k = 0; k < frame_no - corr_time; k++)
            {
                float vector_start[3], vector_stop[3];
                vector_start[0] = traj[k][atom_1_tmp][0] - traj[k][atom_2_tmp][0]; // x-component
                vector_start[1] = traj[k][atom_1_tmp][1] - traj[k][atom_2_tmp][1]; // y-component
                vector_start[2] = traj[k][atom_1_tmp][2] - traj[k][atom_2_tmp][2]; // z-component
                vector_stop[0] = traj[k+corr_time][atom_1_tmp][0] - traj[k+corr_time][atom_2_tmp][0]; // x-component
                vector_stop[1] = traj[k+corr_time][atom_1_tmp][1] - traj[k+corr_time][atom_2_tmp][1]; // y-component
                vector_stop[2] = traj[k+corr_time][atom_1_tmp][2] - traj[k+corr_time][atom_2_tmp][2]; // z-component
                // Length of both vectors -> used to normalize dot product
                float vector_start_length = sqrt(vector_start[0] * vector_start[0] + vector_start[1] * vector_start[1] + vector_start[2] * vector_start[2]);
                float vector_stop_length = sqrt(vector_stop[0] * vector_stop[0] + vector_stop[1] * vector_stop[1] + vector_stop[2] * vector_stop[2]);

                // Calculate dot product of normalized vectors as oacf_tmp; add to oacf[i][1] divided by total number of sampled oacf-values in order to get average
                float oacf_tmp = (vector_start[0] * vector_stop[0] + vector_start[1] * vector_stop[1] + vector_start[2] * vector_stop[2]) / (vector_start_length * vector_stop_length);
                oacf[i][1] += oacf_tmp / (atom_1_no * (frame_no - corr_time));
            }
        }
    }

    // Print performed calculation into output-string (shows auto complete)
    char atom_1_name[3];
    char atom_2_name[3];
    no_to_element(target_atom_1, atom_1_name);
    no_to_element(target_atom_2, atom_2_name);
    if (output != NULL)
    {
        sprintf(output, "\toacf %s %s %i %f\n", atom_1_name, atom_2_name, resolution, timerange);
    }

    return 0;
}
