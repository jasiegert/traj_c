traj:
	clang -lm -Wall dostuff.c read_trajec.c chemistry.c calc/msd.c calc/rdf.c calc/oacf.c -o dostuff.out

