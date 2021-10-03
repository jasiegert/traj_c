#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "read_trajec.h"
#include "chemistry.h"
#include "calc/msd.h"

#define OUTPUTFILE "traj_c.out"

int docalc(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line);

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("Usage: ./dostuff.c xyz-file pbc-file calc-file\n\n");
        return 1;
    }

    // Get file name, atom number and line number{
    char* name = argv[1];
    char* pbc_dat = argv[2];
    char* calc_dat = argv[3];
    int atom_no = get_atoms(name);
    int frame_no = get_lines(name) / (atom_no + 2);
    printf("frame_no: %i\n", frame_no);

    // Read calc-file and dump its content into calc_dump
    FILE *calc_file = fopen(calc_dat, "r");
    if (calc_file == NULL)
    {
        printf("Calc-file %s couldn't be opened.\n", calc_dat);
        return 1;
    }
    fseek (calc_file, 0, SEEK_END);
    int calc_length = ftell(calc_file);
    fseek (calc_file, 0, SEEK_SET);
    char calc_dump[calc_length + 1];
    fread (calc_dump, 1, calc_length, calc_file);
    calc_dump[calc_length] = '\0';
    fclose(calc_file);

    // Create atom array holding atom labels and trajectory array holding coordinates
    float (*traj)[atom_no][3];
    int atom[atom_no];
    traj = malloc(sizeof(*traj)*frame_no);
    if (traj == NULL)
    {
        printf("Couldn't create traj in memory.\n");
        return 1;
    }

    // Read trajectory contents into traj and atom
    if (readxyz(name, frame_no, atom_no, traj, atom) != 0)
    {
        free(traj);
        return 1;
    }
    printf("Trajectory has been read.\n");

    // Read pbc-file into pbc
    printf("Reading pbc-file...");
    float pbc[3][3];
    if (readpbc(pbc_dat, pbc) != 0)
    {
        free(traj);
        return 1;
    }
    printf("done\n");
    
    // Remove com from trajectory
    printf("Removing com...");
    float (*trajcom)[atom_no][3];
    trajcom = malloc(sizeof(*trajcom)*frame_no);
    if (trajcom == NULL)
    {
        printf("Couldn't create trajcom in memory.\n");
        free(traj);
        return 1;
    }
    removecom(frame_no, atom_no, traj, atom, trajcom);
    printf("done\n");



    // Print trajectory to stdout to check it out
    //printarray(frame_no, atom_no, traj, atom);

    // Save trajectory to out.xyz to check it out
    //    char outname[] = "out.xyz";
    //    writexyz(outname, frame_no, atom_no, traj, atom);
    
    // Open output file, write date/time and trajectory/pbc-file
    FILE *output = fopen(OUTPUTFILE, "a");
    time_t current_time = time(NULL);
    fprintf(output, "%s", ctime(&current_time));
    fprintf(output, "Trajectory: %s\nPBC-file: %s\n", name, pbc_dat);
    fclose(output);

    // Parse input
    char* line = strtok(calc_dump, "\n");
    int calc_line_no = 1;
    while (line != NULL)
    {
        printf("\nLine %i: %s\n", calc_line_no, line);
        docalc(frame_no, atom_no, trajcom, atom, line);
        line = strtok(NULL, "\n");
        calc_line_no++;
    }

    free(traj);
    free(trajcom);

    return 0;
}

int docalc(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line)
{
    // printf("Parsing file: %s\n", line);

    // Take first string as calc_name
    char calc_name[10];
    if (sscanf(line, "%s %*[^\n] ", calc_name) == 0)
    {
        //printf("Following line in calc-file does not start with a calculation name:\n\t%s\n\n", line);
        printf("\tDoes not start with a calculation name.\n");
        return 1;
    }
    // Check whether calc_name matches any predefined calculation name, which currently are:
    // msd
    else if (strcmp(calc_name, "msd") == 0)
    {
        // Set default values for resolution and timerange
        printf("\tMSD calculation\n");
        float timerange = 0.4;
        int resolution = 100;
        char target_atom[3];

        // Read atom type from calc-line and complain if atom type is missing or unrecognized
        if ( sscanf(line, "%*s %s %i %f", target_atom, &resolution, &timerange) >= 1)
        {
            int target_atom_no = element_to_no(target_atom);
            if ( target_atom_no == -1)
            {
                printf("\tAtom type %s not recognized.\n", target_atom);
                return 1;
            }
            // Allocate output array msd and string msd_output, then run calculation to fill them with result
            float (*msd)[2];
            msd = malloc(sizeof(msd) * resolution);
            char* msd_output;
            msd_overall(frame_no, atom_no, traj, atom, target_atom_no, resolution, timerange, msd, &msd_output);

            // Write msd into outputcsv and outputstring into outputfile
            char outputcsv[11];
            sprintf(outputcsv, "msd_%s.csv", target_atom);
            savecsv(outputcsv, resolution, 2, msd);
            FILE *output = fopen(OUTPUTFILE, "a");
            fprintf(output, "%s\n", msd_output);
            printf("%s\n", msd_output);

            fclose(output);
            free(msd_output);
            free(msd);
        }
        else
        {
            printf("Following line does not give msd the necessary atom type:\n\t%s\n\n", line);
            return 1;
        }


    }
    else
    {
        printf("\tDoes not start with a valid calculation.\n");
        return 1;
    }


    return 0;
}
