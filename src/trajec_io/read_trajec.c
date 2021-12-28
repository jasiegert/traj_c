#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chemistry.h"
#include "read_trajec.h"

// Reads atom number from file; returns -1 if it fails, otherwise returns atom number
int get_atom_and_frame_no(char *name, int *atom_no, int *frame_no)
{
    FILE *xyz;
    if ( ( xyz = fopen(name, "r")) == NULL)
    {
        printf("xyz-file %s could not be opened!\n", name);
        return 1;
    }

    // Read atom_no from first line of the xyz
    if (fscanf(xyz, " %d", atom_no) != 1)
    {
        printf("Atom number could not be read from the first line of trajectory file!\n");
        fclose(xyz);
        return 1;
    }
    // Count line_no
    rewind(xyz);
    int line_no = 0;
    while (fscanf(xyz, " %*[^\n] \n") != EOF)
    {
        line_no++;
    }

    // Calculate frame_no
    *frame_no = line_no / (*atom_no + 2);
    if (*frame_no * (*atom_no + 2) != line_no)
    {
        printf("xyz-file is not correctly formatted (line_no is not divisible by atom_no + 2).");
        fclose(xyz);
        return 1;
    }

    printf("Found %i atoms and %i frames.\n", *atom_no, *frame_no);
    fclose(xyz);
    return 0;
}

void skipline(FILE *f)
{
    char c;
    do
    {
        c = fgetc(f);
    } while (c != '\n');
}

// Reads trajectoy stored in $name with $frame_no frames and $atom_no atoms; writes coordinates into $traj and atomic numbers into $atom
int readxyz(char *name, int *frame_no_pointer, int *atom_no_pointer, float** trajectory_pointer, int **atom_pointer)
{
    get_atom_and_frame_no(name, atom_no_pointer, frame_no_pointer);
    int atom_no = *atom_no_pointer, frame_no = *frame_no_pointer;

    // Allocate trajectory and atom arrays
    int *atom = malloc(sizeof(int) * atom_no);
    float (*traj)[atom_no][3] = malloc( sizeof(float) * frame_no * atom_no * 3 );
    
    FILE *xyz = fopen(name, "r");
    if (xyz == NULL)
    {
        printf("xyz-file %s could not be read!\n", name);
        return 1;
    }


    // Read atom types from first frame
    char atom_label[3];
    skipline(xyz);
    skipline(xyz);
    for (int j = 0; j < atom_no; j++)
    {
        if (fscanf(xyz, " %s %*[^\n] \n", atom_label) != 1)
        {
            printf("Label of atom %d could not be read!\n", j);
            free(atom);
            free(traj);
            return 1;
        }
        atom[j] = element_to_no(atom_label);
    }
    rewind(xyz);
    printf("Atoms read.\n");

    // Read coordinates into array traj
    for (int i = 0; i < frame_no; i++)
    {
        skipline(xyz);
        skipline(xyz);
        for (int j = 0; j < atom_no; j++)
        {
            if (fscanf(xyz, " %*s %f %f %f \n", &traj[i][j][0], &traj[i][j][1], &traj[i][j][2]) != 3)
            {
                printf("Error reading xyz-file %s in frame %d atom %d!\n", name, i, j);
                free(atom);
                free(traj);
                return 1;
            }
        }
    }
    fclose(xyz);
    printf("Coordinates read.\n");
    *trajectory_pointer = (float*) traj;
    *atom_pointer = atom;
    return 0;
}

int readtraj(char *name, int *frame_no_pointer, int *atom_no_pointer, float **trajectory_pointer, int **atom_pointer)
{
    // Generate name and file pointer for the dat-file
    int namelen = strlen(name);
    char datname[namelen + 1];
    strcpy(datname, name);
    datname[namelen-3] = 'd';
    datname[namelen-2] = 'a';
    datname[namelen-1] = 't';
    
    if ( readdat(datname, frame_no_pointer, atom_no_pointer, trajectory_pointer, atom_pointer) != 0 )
    {
        // Read atom numbers and coordinates from xyz-file
        if (readxyz(name, frame_no_pointer, atom_no_pointer, trajectory_pointer, atom_pointer) != 0)
        {
            printf("Trajectory could not be read from xyz-file.\n");
            return 1;
        }
        int atom_no = *atom_no_pointer, frame_no = *frame_no_pointer;
        printf("Trajectory has been read from xyz-file.\n");
        
        // Write dat-file for future use
        FILE *trajdat;
        if ( (trajdat = fopen(datname, "w")) == NULL)
        {
            printf("Couldn't open %s to write auxiliary dat-file.\n", datname);
        }
        else
        {
            fwrite(atom_no_pointer, 1, sizeof(int), trajdat);
            fwrite(frame_no_pointer, 1, sizeof(int), trajdat);
            fwrite(*atom_pointer, atom_no, sizeof(int), trajdat);
            fwrite(*trajectory_pointer, frame_no * atom_no * 3, sizeof(float), trajdat);
            fclose(trajdat);
            printf(" and written to a dat-file.\n");
        }
    }
 
    return 0;
}

int readdat(char *datname, int *frame_no_pointer, int *atom_no_pointer, float **trajectory_pointer, int **atom_pointer)
{
    FILE* trajdat = fopen(datname, "r");
    if ( trajdat  == NULL )
    {
        printf("Did not find dat-file.\n");
        return 1;
    }
    
    if ( fread(atom_no_pointer, sizeof(int), 1, trajdat) != 1 || fread(frame_no_pointer, sizeof(int), 1, trajdat) != 1)
    {
        printf("Invalid dat-file. Will be overwritten.\n");
        return 1;
    }
    // Allocate trajectory and atom arrays
    int atom_no = *atom_no_pointer, frame_no = *frame_no_pointer;
    *atom_pointer = malloc(sizeof(int) * atom_no);
    *trajectory_pointer = malloc(sizeof(float) * frame_no * atom_no * 3);
    // Read atom numbers and coordinates from dat-file
    if ( fread(*atom_pointer, sizeof(int), atom_no, trajdat) != atom_no || fread(*trajectory_pointer, sizeof(float), frame_no * atom_no * 3, trajdat) != frame_no * atom_no * 3)
    {
        printf("Invalid dat-file. Will be overwritten.\n");
        return 1;
    }

    fclose(trajdat);
    printf("Trajectory has been read from dat-file.\n");
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

