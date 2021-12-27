#include <stdio.h>
#include <stdlib.h>

#include "chemistry.h"

// Reads atom number from file; returns -1 if it fails, otherwise returns atom number
int get_atoms(char *name)
{
    // Open xyz file
    FILE *xyz = fopen(name, "r");
    if (xyz == NULL)
    {
        printf("xyz-file %s could not be read!\n", name);
        return -1;
    }
    // Read atom_no from first line of the xyz
    int atom_no;
    if (fscanf(xyz, " %d", &atom_no) != 1)
    {
        printf("Atom number could not be read from the first line of %s!\n", name);
        fclose(xyz);
        return -1;
    }
    printf("atom_no: %i\n", atom_no);
    fclose(xyz);
    return atom_no;
}

// Count number of lines in file
int get_lines(char *name)
{
    FILE *xyz = fopen(name, "r");
    if (xyz == NULL)
    {
        printf("xyz-file %s could not be read!\n", name);
        return -1;
    }
    // Count lines in file
    int line_no = 0;
    while (fscanf(xyz, " %*[^\n] \n") != EOF)
    {
        line_no++;
    }
    fclose(xyz);
    //printf("line_no: %i\n", line_no);
    return line_no;
}

// Reads trajectoy stored in $name with $frame_no frames and $atom_no atoms; writes coordinates into $traj and atomic numbers into $atom
int readxyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no])
{
    FILE *xyz = fopen(name, "r");
    if (xyz == NULL)
    {
        printf("xyz-file %s could not be read!\n", name);
        return 1;
    }

    // Read atom types from first frame
    char atom_label[3];
    fscanf(xyz, " %*[^\n] \n");
    fscanf(xyz, " %*[^\n] \n");
    for (int j = 0; j < atom_no; j++)
    {
        if (fscanf(xyz, " %s %*[^\n] \n", atom_label) != 1)
        {
            printf("Label of atom %d could not be read!\n", j);
            return 1;
        }
        atom[j] = element_to_no(atom_label);
    }
    rewind(xyz);
    printf("Atoms read.\n");

    // Read coordinates into array traj
    for (int i = 0; i < frame_no; i++)
    {
        fscanf(xyz, " %*[^\n] \n");
        fscanf(xyz, " %*[^\n] \n");
        for (int j = 0; j < atom_no; j++)
        {
            if (fscanf(xyz, " %*s %f %f %f \n", &traj[i][j][0], &traj[i][j][1], &traj[i][j][2]) != 3)
            {
                printf("Error reading xyz-file %s in line %d atom %d!\n", name, i, j);
                return 1;
            }
        }
    }
    fclose(xyz);
    printf("Coordinates read.\n");
    return 0;
}

int readpbc(char *name, float pbc[3][3])
{
    FILE *pbc_dat = fopen(name, "r");
    if (pbc_dat == NULL)
    {
        printf("pbc-file %s could not be read!\n", name);
        return 1;
    }

    for (int i = 0; i < 3; i++)
    {
        if (fscanf(pbc_dat, "%f %f %f \n", &pbc[i][0], &pbc[i][1], &pbc[i][2]) != 3)
        {
            printf("pbc-file %s did not contain 3 readable float values in line %i.\n", name, i);
            return 1;
        }
    }
    return 0;
}

int removecom(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], float trajcom[frame_no][atom_no][3])
{
    // Generate atom masses
    float atom_mass[atom_no];
    float total_mass = 0.0;
    for (int j = 0; j < atom_no; j++)
    {
        atom_mass[j] = no_to_mass(atom[j]);
        total_mass += atom_mass[j];
    }

    // Calculate com for each timestep and subtract it from each elements coord
    for (int i = 0; i < frame_no; i++)
    {
        float com[3] = {0, 0, 0};
        for (int j = 0; j < atom_no; j++)
        {
            com[0] += traj[i][j][0] * atom_mass[j] / total_mass;
            com[1] += traj[i][j][1] * atom_mass[j] / total_mass;
            com[2] += traj[i][j][2] * atom_mass[j] / total_mass;
        }
        for (int j = 0; j < atom_no; j++)
        {
            trajcom[i][j][0] = traj[i][j][0] - com[0];
            trajcom[i][j][1] = traj[i][j][1] - com[1];
            trajcom[i][j][2] = traj[i][j][2] - com[2];
        }
    }

    return 0;
}

// Writes atom labels in $atom and coordinates in $traj as trajectory with $frame_no frames and $atom_no atoms into file $name
int writexyz(char *name, int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no])
{
    // Open output file
    FILE *xyzout = fopen(name, "w+");
    if (xyzout == NULL)
    {
        printf("Couldn't open output xyz-file %s.\n", name);
        return 1;
    }

    // Generate atom labels from atom numbers
    char atom_labels[atom_no][3];
    for (int j = 0; j < atom_no; j++)
    {
        no_to_element(atom[j], atom_labels[j]);
    }

    // Write trajectory content
    for (int i = 0; i < frame_no; i++)
    {
        fprintf(xyzout, "     %i\n\n", atom_no);
        for (int j = 0; j < atom_no; j++)
        {
            fprintf(xyzout, "%3s %20.10f %19.10f %19.10f \n", atom_labels[j], traj[i][j][0], traj[i][j][1], traj[i][j][2]);
        }
    }
    
    // Close output file
    fclose(xyzout);
    return 0;
}

