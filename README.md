# Trajectory analyzer #

This is a tool used to analyze a trajectory, the result of a molecular dynamics (MD) simulation. If you're unfamiliar with this field, you might want to check out the [Theoretical Background](docs/THEORY.md).

# CS50 final project #

Author: Johnny Alexander Jimenez Siegert

City: Halle

Country: Germany

# Installation #

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

- pbc file: example/traj/pbc_CDP_224

    Periodic boundary conditions of the CDP simulation (corresponding to 2x2x4 elementary units of CDP).

- input file: example/calc.inp

 Example input file covering every supported calculation.

Once inside the example/ directory, the calculation can be invoked as:

~~~~~~~~~~~~~~
mkdir results; 
cd results;
../../traj_analyzer ../traj/CDP_224_small.xyz ../traj/pbc_CDP_224 ../calc.inp;
~~~~~~~~~~~~~~

These commands have been collected in example/example.sh, once inside the example/ directory you can run `bash example.sh`. The results will be stored in the directory example/results. For reference the same functions have been computed using the free program package [TRAVIS](http://www.travis-analyzer.de/). The csv-files with the same name in results/ and reference/ should be just about the same (within margin of error).

# Calculations #

Currently three types of calculations are supported: mean square displacement (MSD), radial distribution function (RDF) and orientational auto-correlation function (OACF). 

Each function requires either one (MSD) or two (RDF, OACF) obligatory arguments. The remaining arguments are optional and will be assumed to be their default values, if not specified. All calculations with arguments are:

Calulcation | Parameters | Default | Description
------------|------------|---------|------------
msd         | Atom type  |    -    | The atom type, whose MSD shall be calculated. Only obligatory argument.
^           | Resolution |   100   | Number of correlation times to be sampled. Higher resolution might give a more accurate and less jagged result, but will also increase calculation time.
^           | Time depth |   0.3   | Correlation times will be sampled up until this fraction of the total trajectory length. Higher time depth will give more data by sampling a larger range of correlation times, but the additionally sampled range will be less reliable.
msd_fft     | Atom type  |    -    | The atom type, whose MSD shall be calculated. Only obligatory argument.
^           | Time depth |   0.3   | Correlation times will be sampled up until this fraction of the total trajectory length. Higher time depth will give more data by sampling a larger range of correlation times, but the additionally sampled range will be less reliable.
rdf (rdf_inter)        | Atom type 1 |   -    | First atom type in the rdf-pair.
^ | Atom type 2 |   -    | Second atom type in the rdf-pair.
^           | Bins        |   300  | Number of bins, into which the distances will be sorted. More bins might give slightly more accurate results.
^           | Min distances | 0    | Minimum distance to bin.
^           | Max distances | 4    | Maximum distance to bin.
oacf        | Atom type 1 |   -    | First atom type, which will be the start of the rotating vector.
^           | Atom type 2 |   -    | Second atom type, of which the closest neighbors to the atoms of type 1 will be the end of the rotating vector.
^           | Resolution |   100   | Number of correlation times to be sampled. Higher resolution might give a more accurate and less jagged result, but will also increase calculation time.
^           | Time depth |   0.3   | Correlation times will be sampled up until this fraction of the total trajectory length. Higher time depth will give more data by sampling a larger range of correlation times, but the additionally sampled range will be less reliable.

The input file must contain one line for each calculation to be performed. Each line should contain a calculation name followed by its arguments. Empty lines and lines beginning with # are ignored.

The content of the example input file is as follows:

    msd H 100 0.3
    msd_fft H 0.3
    rdf O O 300 0 4
    rdf_inter O O P 300 0 4
    oacf P O 100 0.

# License #

This project comes with a [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause) license. Feel free to include (whole or in parts), modify or distribute it with your own projects, as long as you retain the information specified in the [license file](../../LICENSE).

I have included parts of the [kissFFT](https://github.com/mborgerding/kissfft) project by Mark Borgerding, which comes with a BSD-3-clause as well. For more information see the [license file](../../src/kissFFT/COPYING) in src/kissfft/.
