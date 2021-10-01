#include <stdio.h>
#include <stdlib.h>
#include "read_trajec.h"
#include "chemistry.h"

int main()
{
    int a1 = element_to_no("O");
    printf("%i\n\n", a1);

    // Get file name, atom number and line number{
    char name[] = "test_traj/CDP_224_example.xyz";
    int atom_no = get_atoms(name);
    int frame_no = get_lines(name) / (atom_no + 2);

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
    readxyz(name, frame_no, atom_no, traj, atom);
    printf("Trajectory has been read.\n");

    // Remove com from trajectory
    float (*trajcom)[atom_no][3];
    trajcom = malloc(sizeof(*traj)*frame_no);
    if (trajcom == NULL)
    {
        printf("Couldn't create trajcom in memory.\n");
        return 1;
    }
    removecom(frame_no, atom_no, traj, atom, trajcom);
    writexyz("out_com.xyz", frame_no, atom_no, trajcom, atom);


    // Print trajectory to stdout to check it out
    //printarray(frame_no, atom_no, traj, atom);
    char outname[] = "out.xyz";
    writexyz(outname, frame_no, atom_no, traj, atom);
    free(traj);
    free(trajcom);
    return 0;
}
