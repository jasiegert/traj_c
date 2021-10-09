#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../chemistry.h"

// Calculates the RDF between two atom types passed as atom_1 and atom_2 in the trajectory traj. Returns RDF by writing into rdf.
int rdf_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], int atom_1, int atom_2, float min_distance, float max_distance, int bin, float rdf[bin][2], char *output)
{
    // Initialize histogram counting all binned atom pair distances
    int histogram[bin];
    for (int i = 0; i < bin; i++)
    {
        histogram[i] = 0;
    }

    // Find all atoms of types atom_1 and atom_2 and save their indices and total number
    int atom_1_index[atom_no], atom_2_index[atom_no];
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

    // Print total number of found atoms of each type
    //char atom_1_name[3], atom_2_name[3];
    //no_to_element(atom_1, atom_1_name);
    //no_to_element(atom_2, atom_2_name);
    //printf("\tFound %i %s atoms and %i %s atoms.\n", atom_1_no, atom_1_name, atom_2_no, atom_2_name);

    // Iterate through each timestep, each atom of type atom_1 and atom of type atom_2
    for (int i = 0; i < frame_no; i++)
    {
        for (int j = 0; j < atom_1_no; j++)
        {
            int atom_1_tmp = atom_1_index[j];
            for (int k = 0; k < atom_2_no; k++)
            {
                int atom_2_tmp = atom_2_index[k];
                // Calculate distance with pbc_dist and if it is between min_distance and max_distance, then add to the appropiate bin
                float dist = pbc_dist(traj[i][atom_1_tmp], traj[i][atom_2_tmp], pbc);
                if ( (dist > min_distance) && (dist < max_distance))
                {
                    int index = floor((dist - min_distance) / (max_distance - min_distance) * bin);
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
        rdf[i][0] = ((i + 0.5) / bin * (max_distance - min_distance) + min_distance) * 100; // average distance for each bin in pm
        rdf[i][1] = (float)histogram[i] / frame_no / ideal_density / shell_volume;  // normalized RDF
    }

    if (output != NULL)
    {
        char atom_1_name[3], atom_2_name[3];
        no_to_element(atom_1, atom_1_name);
        no_to_element(atom_2, atom_2_name);
        sprintf(output, "\trdf %s %s %i %f %f\n\tFound %i %s atoms and %i %s atoms.\n", atom_1_name, atom_2_name, bin, min_distance, max_distance, atom_1_no, atom_1_name, atom_2_no, atom_2_name);
    }

    return 0;
}

// Calculates the RDF between two atom types passed as atom_1 and atom_2 in the trajectory traj. Returns RDF by writing into rdf.
int rdf_intermolecular(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], int atom_1, int atom_2, int atom_central, float min_distance, float max_distance, int bin, float rdf[bin][2], char *output)
{
    // Allow calculation only if atom types are H, O or the central atom
    if ( !(((atom_1 == 1) || (atom_1 == 8) || (atom_1 == atom_central)) && ((atom_2 == 1) || (atom_2 == 8) || (atom_2 == atom_central))))
    {
        printf("\tRDF with central atom is only supported for atom types H, O or the central atom.");
        return 1;
    }

    // Initialize histogram counting all binned atom pair distances
    int histogram[bin];
    for (int i = 0; i < bin; i++)
    {
        histogram[i] = 0;
    }

    // Build molecule index to identify atoms of the same molecule
    int molecule[atom_no];
    // Number central atoms in order of their appearance in atom / trajectory
    int atom_central_index[atom_no], atom_O_index[atom_no], atom_H_index[atom_no];
    int atom_central_no = 0, atom_O_no = 0, atom_H_no = 0;
    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_central)
        {
            molecule[i] = atom_central_no;
            atom_central_index[atom_central_no] = i;
            atom_central_no++;
        }
    }
    // Assign oxygen to the same molecule as their closest central atom
    int atom_O = element_to_no("O");
    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_O)
        {
            int neighbor_central = nextneighbor_in_traj(traj[0][i], atom_no, traj[0], atom, atom_central, pbc);
            atom_O_index[atom_O_no] = i;
            atom_O_no++;
            molecule[i] = molecule[neighbor_central];
        }
    }
    // Assign hydrogen to the same molecule as their closest oxygen
    int atom_H = element_to_no("H");
    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_H)
        {
            int neighbor_O = nextneighbor_in_traj(traj[0][i], atom_no, traj[0], atom, atom_O, pbc);
            atom_H_index[atom_H_no] = i;
            atom_H_no++;
            molecule[i] = molecule[neighbor_O];
        }
    }

    int atom_1_no, atom_2_no;
    int atom_1_index[atom_no], atom_2_index[atom_no];
    if (atom_1 == 1)
        {memcpy(atom_1_index, atom_H_index, sizeof(atom_1_index)); atom_1_no = atom_H_no;}
    else if (atom_1 == 8)
        {memcpy(atom_1_index, atom_O_index, sizeof(atom_1_index)); atom_1_no = atom_O_no;}
    else if (atom_1 == atom_central)
        {memcpy(atom_1_index, atom_central_index, sizeof(atom_1_index)); atom_1_no = atom_central_no;}
    else
        return 1;
    if (atom_2 == 1)
        {memcpy(atom_2_index, atom_H_index, sizeof(atom_2_index)); atom_2_no = atom_H_no;}
    else if (atom_2 == 8)
        {memcpy(atom_2_index, atom_O_index, sizeof(atom_2_index)); atom_2_no = atom_O_no;}
    else if (atom_2 == atom_central)
        {memcpy(atom_2_index, atom_central_index, sizeof(atom_2_index)); atom_2_no = atom_central_no;}
    else
        return 1;
    
    // Print total number of found atoms of each type
    //char atom_1_name[3], atom_2_name[3], atom_central_name[3];
    //no_to_element(atom_1, atom_1_name);
    //no_to_element(atom_2, atom_2_name);
    //no_to_element(atom_central, atom_central_name);
    //printf("\tFound %i %s atoms and %i %s atoms.\n", atom_1_no, atom_1_name, atom_2_no, atom_2_name);

    // Iterate through each timestep, each atom of type atom_1 and atom of type atom_2
    for (int i = 0; i < frame_no; i++)
    {
        for (int j = 0; j < atom_1_no; j++)
        {
            int atom_1_tmp = atom_1_index[j];
            for (int k = 0; k < atom_2_no; k++)
            {
                int atom_2_tmp = atom_2_index[k];
                if (molecule[atom_2_tmp] != molecule[atom_1_tmp])
                {
                    // Calculate distance with pbc_dist and if it is between min_distance and max_distance, then add to the appropiate bin
                    float dist = pbc_dist(traj[i][atom_1_tmp], traj[i][atom_2_tmp], pbc);
                    if ( (dist > min_distance) && (dist < max_distance))
                    {
                        int index = floor((dist - min_distance) / (max_distance - min_distance) * bin);
                        histogram[index] += 1;
                    }
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
        rdf[i][0] = ((i + 0.5) / bin * (max_distance - min_distance) + min_distance) * 100; // average distance for each bin in pm
        rdf[i][1] = (float)histogram[i] / frame_no / ideal_density / shell_volume;  // normalized RDF
    }

    // Print used parameters into output-string
    if (output != NULL)
    {
        char atom_1_name[3], atom_2_name[3], atom_central_name[3];
        no_to_element(atom_1, atom_1_name);
        no_to_element(atom_2, atom_2_name);
        no_to_element(atom_central, atom_central_name);
        sprintf(output, "\trdf %s %s %i %f %f\n\tFound %i %s atoms, %i %s atoms and %i %s atoms.\n", atom_1_name, atom_2_name, bin, min_distance, max_distance, atom_central_no, atom_central_name, atom_1_no, atom_1_name, atom_2_no, atom_2_name);
    }

    return 0;
}
