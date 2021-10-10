#ifndef MSD_H
#define MSD_H

int msd_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, int resolution, float timerange, float msd[resolution][2], char *output);
int msd_fft(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, float timerange, float msd[frame_no][2], char *output);

#endif /* MSD_H */
