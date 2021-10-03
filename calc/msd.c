#include <math.h>
#include <stdio.h>

int msd_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, int resolution, float timerange, float msd[resolution][2], char **output)
{
    float coord_diff;
    int sampled;
    float timestep = 0.5;

    // Count number of target atoms in trajectory
    int target_atom_no = 0;
    for (int j = 0; j < atom_no; j++)
    {
        if (atom[j] == target_atom)
        {
            target_atom_no++;
        }
    }

    // Create correlation time for each resolution point 
    for (int i = 1; i <= resolution; i++)
    {
        int corr_time = round(i / (float)resolution * timerange * frame_no);
        msd[i][0] = corr_time * timestep / 1000;   // Correlation time in ps
        msd[i][1] = 0;             // MSD for correlation time in pm^2

        // Check each atom and calculate msd if it is of target type
        for (int j = 0; j < atom_no; j++)
        {
            if (atom[j] == target_atom)
            {
                // Iterate through each correlation window (which are just shifted by one timestep)
                for (int k = 0; k < frame_no - corr_time; k++)
                {
                    coord_diff = (traj[k + corr_time][j][0] - traj[k][j][0]) * (traj[k + corr_time][j][0] - traj[k][j][0]);     // x component
                    coord_diff += (traj[k + corr_time][j][1] - traj[k][j][1]) * (traj[k + corr_time][j][1] - traj[k][j][1]);    // y component
                    coord_diff += (traj[k + corr_time][j][2] - traj[k][j][2]) * (traj[k + corr_time][j][2] - traj[k][j][2]);    // z component
                    sampled = target_atom_no * (frame_no - corr_time);
                    coord_diff /= sampled;
                    msd[i][1] += coord_diff * 1E4; //pm^2
                }
            }
        }
    }

    return 0;
}
