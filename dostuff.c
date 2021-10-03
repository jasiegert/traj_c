#include <stdio.h>
#include <stdlib.h>
#include "read_trajec.h"
#include "chemistry.h"
#include "calc/msd.h"

int docalc(char* calcline);

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
    printf("Removing com...\n");
    float (*trajcom)[atom_no][3];
    trajcom = malloc(sizeof(*trajcom)*frame_no);
    if (trajcom == NULL)
    {
        printf("Couldn't create trajcom in memory.\n");
        free(traj);
        return 1;
    }
    printf("running removecom\n");
    removecom(frame_no, atom_no, traj, atom, trajcom);
    printf("done\n");

//    printf("Writing to out.xyz...");
//    writexyz("out_com.xyz", frame_no, atom_no, trajcom, atom);
//    printf("done\n");



    // Print trajectory to stdout to check it out
    //printarray(frame_no, atom_no, traj, atom);

    // Save trajectory to out.xyz to check it out
//    char outname[] = "out.xyz";
//    writexyz(outname, frame_no, atom_no, traj, atom);
//    free(traj);
//    free(trajcom);

    // Calculate msd
    printf("Calculationg msd...");
    int resolution = 100;
    float (*msd)[2];
    msd = malloc(sizeof(msd) * resolution);
//    float msd[resolution][2];
    char* msd_output;
    msd_overall(frame_no, atom_no, trajcom, atom, 1, resolution, 0.4, msd, &msd_output);
    printf("done\nWriting to csv...");
    savecsv("msd_H.csv", resolution, 2, msd);
    printf("done\n");

    return 0;
}
