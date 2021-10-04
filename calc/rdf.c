#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../chemistry.h"

int rdf_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], int atom_1, int atom_2, float min_distance, float max_distance, int bin, float rdf[bin][2], char *output)
{
    int histogram[bin];
    for (int i = 0; i < bin; i++)
    {
        histogram[i] = 0;
    }

    int atom_1_index[atom_no];
    int atom_2_index[atom_no];
    int atom_1_no = 0, atom_2_no = 0;


    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_1)
        {
            atom_1_index[atom_1_no] = i;
            atom_1_no++;
        }

        if (atom[i] == atom_2)
        {
            atom_2_index[atom_2_no] = i;
            atom_2_no++;
        }
    }
    atom_1_index[atom_1_no] = -1;
    atom_2_index[atom_2_no] = -1;
    char atom_1_name[3], atom_2_name[3];
    no_to_element(atom_1, atom_1_name);
    no_to_element(atom_2, atom_2_name);
    printf("\tFound %i %s atoms and %i %s atoms.\n", atom_1_no, atom_1_name, atom_2_no, atom_2_name);

    for (int i = 0; i < frame_no; i++)
    {
        for (int j = 0; j < atom_1_no; j++)
        {
            int atom_1_tmp = atom_1_index[j];
            for (int k = 0; k < atom_2_no; k++)
            {
                int atom_2_tmp = atom_2_index[k];
                float dist = pbc_dist(traj[i][atom_1_tmp], traj[i][atom_2_tmp], pbc);
                if ( (dist > min_distance) && (dist < max_distance))
                {
                    int index = floor((dist - min_distance) / (max_distance - min_distance) * bin);
                    //rdf[index][1] = rdf[index][1] + 1;
                    histogram[index] += 1;
                }
            }
        }
    }
    

    // Calculate ideal density and shell volumes in order to normalize RDF
     float pbc_volume = pbc[0][0] * pbc[1][1] * pbc[2][2]; // only valid for orthogonal cell -> otherwise triple-dot product (a x b) * c
     float ideal_density = atom_1_no * atom_2_no / pbc_volume;
    for (int i = 0; i < bin; i++)
    {
        float r_min = (i / (float)bin * (max_distance - min_distance) + min_distance);
        float r_max = ((i + 1) / (float)bin * (max_distance - min_distance) + min_distance);
        float shell_volume = 4 * M_PI / 3 * (r_max * r_max * r_max - r_min * r_min * r_min);
        rdf[i][0] = ((i + 0.5) / bin * (max_distance - min_distance) + min_distance) * 100;
        rdf[i][1] = (float)histogram[i] / frame_no / ideal_density / shell_volume;
    }



    if (output != NULL)
    {
        sprintf(output, "");
    }


    return 0;
}
