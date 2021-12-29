# Trajectory analyzer #

This is a tool used to analyze a trajectory, the result of a molecular dynamics (MD) simulation.

# Installtion #

In order to build the executable, simply invoke `make`. This will build all the necessary object files and an executable traj_analyzer. You can remove the object files with `make clean`.
Compilation can be sped up using the `-j N` option, where N is the number of CPU cores. On a system with 4 cores, you could therefore run:

    make -j 4
    make clean

If you wish to remove all compiled files including the executable, run `make remove`.

# Usage #

The executable traj_analyzer expects three arguments:
- Trajectory (.xyz)
 
 A trajectory, usually produced by a molecular dynamics software. One timestep is assumed to represent 0.5 fs.  Currently only the xyz file format is supported.

- pbc file

 A text file containing the periodic boundary conditions. The first three lines must contain three numbers seperated by a space each.

- input file

 A text file containing on each line the name of a calculation followed by its parameters.

With these three, you can invoke the program as:

    traj_analyzer trajectory.xyz pbc.txt input.txt

# Example #

A small example is included in the directory example/. It contains the necessary inputs for a simple calculation:
- Trajectory: example/traj/CDP_224.xyz

 Trajectory of cesium dihydrogen phosphate (CDP) over 1000 timesteps. 
