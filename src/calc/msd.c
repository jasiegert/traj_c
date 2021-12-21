#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "../chemistry.h"
#include "../kissFFT/kiss_fftr.h"

int msd_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], int atom[atom_no], int target_atom, int resolution, float timerange, float msd[resolution][2], char *output)
{
    float coord_diff;
    int sampled;
    float timestep = 0.5; // in fs

    // Count number of target atoms in trajectory
    int target_atom_no = 0;
    for (int j = 0; j < atom_no; j++)
    {
        if (atom[j] == target_atom)
        {
            target_atom_no++;
        }
    }

    // Iterate through sampled (resolution + 1) correlation times (including 0)
    for (int i = 0; i < resolution + 1; i++)
    {
        // Create correlation time for each resolution point 
        int corr_time = round(i / (float)resolution * timerange * frame_no);
        msd[i][0] = corr_time * timestep / 1000;   // Correlation time in ps
        msd[i][1] = 0;             // MSD for correlation time in pm^2

        // Check each atom and calculate msd if it is of target type
        for (int j = 0; j < atom_no; j++)
        {
            if (atom[j] == target_atom)
            {
                // Iterate through each correlation window (which are just shifted by one timestep)
                for (int k = 0; k < frame_no - corr_time; k++)
                {
                    coord_diff = (traj[k + corr_time][j][0] - traj[k][j][0]) * (traj[k + corr_time][j][0] - traj[k][j][0]);     // x component
                    coord_diff += (traj[k + corr_time][j][1] - traj[k][j][1]) * (traj[k + corr_time][j][1] - traj[k][j][1]);    // y component
                    coord_diff += (traj[k + corr_time][j][2] - traj[k][j][2]) * (traj[k + corr_time][j][2] - traj[k][j][2]);    // z component
                    sampled = target_atom_no * (frame_no - corr_time);
                    coord_diff /= sampled;
                    msd[i][1] += coord_diff * 1E4; //pm^2
                }
            }
        }
    }

    // Apply linear fit to find D
    float slope, intercept, R;
    linregress_array(resolution + 1, msd, 0.5, 1, &slope, &intercept, &R);
    float D = slope / 6;

    char atom_name[3];
    no_to_element(target_atom, atom_name);
    if (output != NULL)
    {
        sprintf(output, "\tmsd %s %i %f\n\tFound %i %s atoms.\n\tDiffusion coefficient: %f (R^2 = %f)\n", atom_name, resolution, timerange, target_atom_no, atom_name,  D, R * R);
    }

    return 0;
}

// Computes MSD using FFT, for the algorithm see:
//https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft
int msd_fft(int frame_no, int atom_no, const float traj[frame_no][atom_no][3], const int atom[atom_no], int target_atom, float timerange, float msd[][2], char *output)
{
    // Alloc buffers for the real FFT and real inverse FFT
    kiss_fftr_cfg fft_buffer = kiss_fftr_alloc(2 * frame_no,0,NULL,NULL);
    kiss_fftr_cfg fft_inv_buffer = kiss_fftr_alloc(2 * frame_no,1,NULL,NULL);
    
    // Initialize msd array for all correlation times to 0
    int msd_len = round(timerange * frame_no) + 1;
    for (int i = 0; i < msd_len; i++)
    {
        msd[i][0] = i * 0.5 / 1E3; // correlation time in ps
        msd[i][1] = 0;
    }

    // Find number of target_atoms
    int target_atom_no = 0;
    for (int j = 0; j < atom_no; j++)
    {
        if (atom[j] == target_atom)
        {
            target_atom_no++;
        }
    }

    // Check each atom whether or not it is of target_type
    for (int j = 0; j < atom_no; j++)
    {
        if (atom[j] == target_atom)
        {
            // coords_for_fft will hold zero-padded coordinates and ultimately S2 (result of inverse FFT)
            kiss_fft_scalar coords_for_fft[3][frame_no * 2];
            // fft_results will hold fourier transform of coordinates
            kiss_fft_cpx fft_results[3][frame_no + 1];
            
            // Write coordinates of atom j into coords_for_fft
            for (int i = 0; i < frame_no; i++)
            {
               coords_for_fft[0][i] = traj[i][j][0];
               coords_for_fft[1][i] = traj[i][j][1];
               coords_for_fft[2][i] = traj[i][j][2];
            }
            // Zero-pad data to get acyclic correlation
            for (int i = frame_no; i < frame_no * 2; i++)
            {
               coords_for_fft[0][i] = 0;
               coords_for_fft[1][i] = 0;
               coords_for_fft[2][i] = 0;
            }

            // Compute S2 (autocorrelation, see stackexchange)
            float S2[frame_no];
            for (int k = 0; k < 3; k++)
            {
                // Fourier transformation into fft_results
                kiss_fftr(fft_buffer, coords_for_fft[k], fft_results[k]);
                // Element-wise fft_results * fft_results.conjugate(); scale with 1/nfft to reverse scaling by kissFFT
                for (int i = 0; i < frame_no + 1; i++)
                {
                    fft_results[k][i].r = fft_results[k][i].r * fft_results[k][i].r + fft_results[k][i].i * fft_results[k][i].i;
                    fft_results[k][i].r *= 1.0 / (2 * frame_no);
                    fft_results[k][i].i = 0;
                }
                // Inverse Fourier transformation to get correlation
                kiss_fftri(fft_inv_buffer, fft_results[k], coords_for_fft[k]);
            }
            for (int i = 0; i < frame_no; i++)
            {
                // Sum over three dimensions scaled with frame_no - corr_time = number of sampled time frames
                S2[i] = (coords_for_fft[0][i] + coords_for_fft[1][i] + coords_for_fft[2][i]) / (frame_no - i);
            }

            // Compute S1 (see stackexchange)
            float S1[frame_no];
            float D[frame_no];
            float Q = 0;
            for (int i = 0; i < frame_no; i++)
            {
                D[i] = traj[i][j][0] * traj[i][j][0] + traj[i][j][1] * traj[i][j][1] + traj[i][j][2] * traj[i][j][2];
                Q += 2 * D[i];
            }
            S1[0] = Q / frame_no;
            for (int i = 1; i < frame_no; i++)
            {
                Q -= D[i-1] + D[frame_no - i];
                S1[i] = Q / (frame_no - i);
            }

            // Add msd of the atom (divided by atom no to get average) to msd-array
            for (int i = 0; i < msd_len; i++)
            {
                msd[i][1] += (S1[i] - 2 * S2[i]) / target_atom_no * 1E4; // msd in pm^2
            }
        }
    }

    // Apply linear regression to find D
    float slope, intercept, R;
    linregress_array(msd_len, msd, 0.5, 1, &slope, &intercept, &R);
    float D = slope / 6;

    // Put output string into *output containing passed parameters, number of atoms and diffusion coefficient
    char atom_name[3];
    no_to_element(target_atom, atom_name);
    if (output != NULL)
    {
        sprintf(output, "\tmsd_fft %s %f\n\tFound %i %s atoms.\n\tDiffusion coefficient: %f (R^2 = %f)\n", atom_name, timerange, target_atom_no, atom_name,  D, R * R);
    }

    free(fft_buffer);
    free(fft_inv_buffer);
    return 0;
}
