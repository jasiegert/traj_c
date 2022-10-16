/**
* @file dostuff.h
* @author Johnny Alexander Jimenez Siegert
* @brief Main control unit of the project, which is compiled into the executable.
*/


/**
* @brief Parses a line from the input file and dispatches appropriate calculation.
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj coordinates of the trajectory as a 3D array.
* @param[in] pbc periodic boundary conditions with one cell vector in each row.
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] line input line to be parsed.
* @return int -1 if the line is to be ignored (empty or comment), 1 if a calculation fails and otherwise 0
*/
int docalc(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line
            );

/**
* @brief Handles MSD calculation.
* @param[in] frame_no number of frames in trajectory.
* @param[in] atom_no number of atoms in trajectory.
* @param[in] traj coordinates of the trajectory as a 3D array.
* @param[in] atom atomic numbers of all atoms in the trajectory.
* @param[in] line input line to be parsed.
*
* Periodic boundary conditions (pbc) are not used, but passed anyway to fit in with the other calc_* functions.
*
*/
int calc_msd(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line
            );
/**
* @brief Handles MSD calculation via fast-fourier transformation.
* @copydetails calc_msd()
*
* Periodic boundary conditions (pbc) are not used, but passed anyway to fit in with the other calc_* functions.
*
*/
int calc_msd_fft(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], char* line
            );
/**
* @brief Handles RDF calculation.
* @copydetails docalc()
*/
int calc_rdf(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line
            );
/**
* @brief Handles intermolecular RDF calculation.
* @copydetails docalc()
*/
int calc_rdf_inter(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line
            );
/**
* @brief Handles OACF calculation.
* @copydetails docalc()
*/
int calc_oacf(
            int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], char* line
            );

/**
* @brief Appends the passed string to the outputfile and prints it to stdout.
*
* The outputfile is defined in the implementation via the Macro OUTPUTFILE.
*
*/
void appendoutput(char outputstring[]);
