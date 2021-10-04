#ifndef RDF_H
#define RDF_H

int rdf_overall(int frame_no, int atom_no, float traj[frame_no][atom_no][3], float pbc[3][3], int atom[atom_no], int atom_1, int atom_2, float min_distance, float max_distance, int bin, float rdf[bin][2], char *output);

#endif /* RDF_H */
