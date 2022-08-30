#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "dostuff.h"
#include "trajec_io/read_trajec.h"
#include "trajec_io/chemistry.h"
#include "calc/mathtools.h"
#include "calc/msd.h"
#include "calc/rdf.h"
#include "calc/oacf.h"

#define OUTPUTFILE "traj_analysis.out"
#define OUTSTRINGLENGTH 100
#define XYZ_NAME_MAX 20


int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("Usage: traj_analyzer xyz-file pbc-file calc-file\n");
        return 1;
    }

    // Get file name, atom number and line number{
    char* name = argv[1];
    char* pbc_dat = argv[2];
    char* calc_dat = argv[3];

    // Read calc-file and dump its content into calc_dump
    FILE *calc_file = fopen(calc_dat, "r");
    if (calc_file == NULL)
    {
        printf("Calc-file %s couldn't be opened.\n", calc_dat);
        return 1;
    }
    fseek (calc_file, 0, SEEK_END);
    size_t calc_length = ftell(calc_file);
    fseek (calc_file, 0, SEEK_SET);
    char calc_dump[calc_length + 1];
    if (fread (calc_dump, 1, calc_length, calc_file) != calc_length)
    {
        printf("Couldn't read calc-file %s!\n", calc_dat);
        return 1;
    }
    calc_dump[calc_length] = '\0';
    fclose(calc_file);

    // Read pbc-file into pbc
    printf("Reading pbc-file...");
    float pbc[3][3];
    if (readpbc(pbc_dat, pbc) != 0)
    {
        printf("Couldn't read pbc-file.\n");
        return 1;
    }
    printf("done\n");

    clock_t c_start, c_end;
    c_start = clock();

    // Read trajectory contents into traj and atom
    int frame_no, atom_no;
    int *atom = NULL;
    float *trajectory_as_pointer = NULL;
    if (readtraj(name, &frame_no, &atom_no, &trajectory_as_pointer, &atom) != 0)
    {
        free(atom);
        free(trajectory_as_pointer);
        return 1;
    }
    float (*traj)[atom_no][3] = ( float (*)[atom_no][3] ) trajectory_as_pointer;
    c_end = clock();
    printf("\ttook %f s in CPU time.\n", (double)(c_end - c_start) / CLOCKS_PER_SEC);

    // Remove com from trajectory
    printf("Removing com...");
    float (*trajcom)[atom_no][3];
    trajcom = malloc(sizeof(*trajcom)*frame_no);
    if (trajcom == NULL)
    {
        printf("Couldn't create trajcom in memory.\n");
        free(traj);
        free(atom);
        return 1;
    }
    removecom(frame_no, atom_no, traj, atom, trajcom);
    printf("done\n");
    
    // Open output file, write date/time and trajectory/pbc-file
    FILE *output = fopen(OUTPUTFILE, "a");
    time_t current_time = time(NULL);
    fprintf(output, "%s", ctime(&current_time));
    fprintf(output, "Trajectory: %s\nPBC-file: %s\n", name, pbc_dat);
    fclose(output);

    // Parse input line by line and perform calculation
    char* line = strtok(calc_dump, "\n");
    while (line != NULL)
    {
        c_start = clock();
        int calc_status = docalc(frame_no, atom_no, trajcom, pbc, atom, line);
        if ( calc_status == 0)
        {
            c_end = clock();
            printf("\ttook %f s in CPU time.\n\n", (double)(c_end - c_start) / CLOCKS_PER_SEC);
        }
        else if (calc_status == -1)
        {
            // do nothing, most likely because it's an empty line or a comment
        }
        else
        {
            printf("\tsomething went wrong.\n\n");
        }
        line = strtok(NULL, "\n");
    }

    free(traj);
    free(trajcom);
    free(atom);
    return 0;
}

int docalc(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line)
{
    // Take first string as calc_name, skip lines that are empty or start with '#'
    char calc_name[10];
    if ( (sscanf(line, "%s %*[^\n] ", calc_name) == 0) || ( calc_name[0] == '#'))
    {
        return -1;
    }
    printf("%s\n", line);

    // Check whether calc_name matches any predefined calculation name, which currently are:
    // msd, msd_fft, rdf, rdf_inter, oacf, removecom
    if (strcmp(calc_name, "msd") == 0)
    {
        return calc_msd(frame_no, atom_no, traj, atom, line);
    }
    else if (strcmp(calc_name, "msd_fft") == 0)
    {
        return calc_msd_fft(frame_no, atom_no, traj, atom, line);
    }
    else if (strcmp(calc_name, "rdf") == 0)
    {
        return calc_rdf(frame_no, atom_no, traj, pbc, atom, line);
    }
    else if (strcmp(calc_name, "rdf_inter") == 0)
    {
        return calc_rdf_inter(frame_no, atom_no, traj, pbc, atom, line);
    }
    else if (strcmp(calc_name, "oacf") == 0)
    {
        return calc_oacf(frame_no, atom_no, traj, pbc, atom, line);
    }
    else if (strcmp(calc_name, "removecom") == 0)
    {
        char output_xyz_name[XYZ_NAME_MAX];
        if ( sscanf(line, "%*s %s", output_xyz_name) >= 1)
        {
            printf("\tSaving trajectory with removed com to %s\n", output_xyz_name);
            writexyz(output_xyz_name, frame_no, atom_no, traj, atom);
            return 0;
        }
        else
        {
            printf("\tPlease provide a name for removecom.\n");
            return 1;
        }
    }
    else
    {
        printf("\tDoes not start with a valid calculation.\n");
        return 1;
    }
}


int calc_msd(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line)
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
        if (target_atom_no == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom);
            return 1;
        }
        
        // Allocate output array msd and string msd_output, then run calculation to fill them with result
        int msd_len = resolution + 1;
        float (*msd)[2];
        msd = malloc( sizeof(msd) * msd_len );
        char msd_output[OUTSTRINGLENGTH];
        msd_overall(frame_no, atom_no, traj, atom, target_atom_no, resolution, timerange, msd, msd_output);

        // Write msd into outputcsv and outputstring into outputfile
        char outputcsv[11];
        sprintf(outputcsv, "msd_%s.csv", target_atom);
        char headerstring[] = "# time / ps\t\t\tMSD / A²";
        savecsv(outputcsv, msd_len, 2, msd, headerstring);

        appendoutput(msd_output);

        free(msd);
        return 0;
    }
    else
    {
        printf("Does not contain the necessary atom type for MSD calculation.\n");
        return 1;
    }
}

int calc_msd_fft(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line)
{
    // Set default values for resolution and timerange
    printf("\tMSD calculation (FFT)\n");
    float timerange = 0.3;
    char target_atom[3];

    // Read atom type from calc-line and complain if atom type is missing or unrecognized
    if ( sscanf(line, "%*s %s %f", target_atom, &timerange) >= 1)
    {
        int target_atom_no = element_to_no(target_atom);
        if ( target_atom_no == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom);
            return 1;
        }
        // Allocate output array msd and string msd_output, then run calculation to fill them with result
        int msd_len = round(frame_no * timerange) + 1;
        float (*msd)[2];
        msd = calloc(msd_len, sizeof(msd));
        char msd_output[OUTSTRINGLENGTH];
        msd_fft(frame_no, atom_no, traj, atom, target_atom_no, timerange, msd, msd_output);

        // Write msd into outputcsv and outputstring into outputfile
        char outputcsv[15];
        sprintf(outputcsv, "msd_fft_%s.csv", target_atom);
        char headerstring[] = "# time / ps\t\t\tMSD / A²";
        savecsv(outputcsv, msd_len, 2, msd, headerstring);
        
        appendoutput(msd_output);

        free(msd);
        return 0;
    }
    else
    {
        printf("Does not contain the necessary atom type for OACF.\n");
        return 1;
    }
}

int calc_rdf(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line)
{
    // Set default values for resolution and timerange
    printf("\tRDF calculation\n");
    int bin = 100;
    float min_distance = 0.0;
    float max_distance = 4.0;
    char target_atom_1[3];
    char target_atom_2[3];

    // Read atom type from calc-line and complain if atom type is missing or unrecognized
    if ( sscanf(line, "%*s %s %s %i %f %f", target_atom_1, target_atom_2, &bin, &min_distance, &max_distance) >= 2)
    {
        int target_atom_no_1 = element_to_no(target_atom_1);
        int target_atom_no_2 = element_to_no(target_atom_2);
        if ( target_atom_no_1 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_1);
            return 1;
        }
        if ( target_atom_no_2 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_2);
            return 1;
        }
        // Allocate output array msd and string msd_output, then run calculation to fill them with result
        float (*rdf)[2];
        rdf = malloc(sizeof(rdf) * bin);
        char rdf_output[OUTSTRINGLENGTH];
        rdf_overall(frame_no, atom_no, traj, pbc, atom, target_atom_no_1, target_atom_no_2, min_distance, max_distance, bin, rdf, rdf_output);

        // Write msd into outputcsv and outputstring into outputfile
        char outputcsv[14];
        sprintf(outputcsv, "rdf_%s_%s.csv", target_atom_1, target_atom_2);
        char headerstring[] = "# d / A\t\t\t\tg(r)";
        savecsv(outputcsv, bin, 2, rdf, headerstring);

        appendoutput(rdf_output);

        free(rdf);
        return 0;
    }
    else
    {
        printf("Does not contain the necessary two atom types for RDF.\n");
        return 1;
    }
}

int calc_rdf_inter(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line)
{
    // Set default values for resolution and timerange
    printf("\tRDF calculation\n");
    int bin = 100;
    float min_distance = 0.0;
    float max_distance = 4.0;
    char target_atom_1[3];
    char target_atom_2[3];
    char central_atom[3];

    // Read atom type from calc-line and complain if atom type is missing or unrecognized
    if ( sscanf(line, "%*s %s %s %s %i %f %f", target_atom_1, target_atom_2, central_atom, &bin, &min_distance, &max_distance) >= 3)
    {
        int target_atom_no_1 = element_to_no(target_atom_1);
        int target_atom_no_2 = element_to_no(target_atom_2);
        int central_atom_no = element_to_no(central_atom);
        if ( target_atom_no_1 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_1);
            return 1;
        }
        if ( target_atom_no_2 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_2);
            return 1;
        }
        // Allocate output array msd and string msd_output, then run calculation to fill them with result
        float (*rdf)[2];
        rdf = malloc(sizeof(rdf) * bin);
        char rdf_output[OUTSTRINGLENGTH];
        rdf_intermolecular(frame_no, atom_no, traj, pbc, atom, target_atom_no_1, target_atom_no_2, central_atom_no, min_distance, max_distance, bin, rdf, rdf_output);

        // Write msd into outputcsv and outputstring into outputfile
        char outputcsv[20];
        sprintf(outputcsv, "rdf_inter_%s_%s.csv", target_atom_1, target_atom_2);
        char headerstring[] = "# d / A\t\t\t\tg(r)";
        savecsv(outputcsv, bin, 2, rdf, headerstring);
        
        appendoutput(rdf_output);

        free(rdf);
        return 0;
    }
    else
    {
        printf("Does not contain the necessary three atom types for intermolecular RDF.\n");
        return 1;
    }
}

int calc_oacf(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line)
{
    // Set default values for resolution and timerange
    printf("\tOACF calculation\n");
    float timerange = 0.4;
    int resolution = 100;
    char target_atom_1[3];
    char target_atom_2[3];

    // Read atom type from calc-line and complain if atom type is missing or unrecognized
    if ( sscanf(line, "%*s %s %s %i %f", target_atom_1, target_atom_2, &resolution, &timerange) >= 2)
    {
        int target_atom_no_1 = element_to_no(target_atom_1);
        int target_atom_no_2 = element_to_no(target_atom_2);
        if ( target_atom_no_1 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_1);
            return 1;
        }
        if ( target_atom_no_2 == -1)
        {
            printf("\tAtom type %s not recognized.\n", target_atom_2);
            return 1;
        }

        // Allocate output array oacf and string oacf_output, then run calculation to fill them with result
        float (*oacf)[2];
        oacf = malloc(sizeof(oacf) * resolution);
        char oacf_output[OUTSTRINGLENGTH];
        oacf_overall(frame_no, atom_no, traj, atom, pbc, target_atom_no_1, target_atom_no_2, resolution, timerange, oacf, oacf_output);

        // Write msd into outputcsv and outputstring into outputfile
        char outputcsv[15];
        sprintf(outputcsv, "oacf_%s_%s.csv", target_atom_1, target_atom_2);
        char headerstring[] = "# time / ps\t\t\tg(r)";
        savecsv(outputcsv, resolution, 2, oacf, headerstring);

        appendoutput(oacf_output);

        free(oacf);
        return 0;
    }
    else
    {
        printf("Does not contain the necessary two atom types for OACF.\n");
        return 1;
    }
}

void appendoutput(char outputstring[])
{
    FILE *output = fopen(OUTPUTFILE, "a");
    if (output == NULL)
    {
        return;
    }
    fprintf(output, "%s\n", outputstring);
    printf("%s", outputstring);
    fclose(output);
}
