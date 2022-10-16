#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "chemistry.h"

//typedef struct atominfo
//{
//    int atom_no;
//    float atom_mass;
//    char symbol[3];
//} atominfo;

int no_to_atominfo(int atom_no, atominfo *this_atom)
{
    switch (atom_no)
    {
        case 1: this_atom->atom_no = 1; this_atom->atom_mass = 1.00798; strcpy(this_atom->symbol, "H");  break;
        case 2: this_atom->atom_no = 2; this_atom->atom_mass = 4.0026; strcpy(this_atom->symbol, "He");  break;
        case 3: this_atom->atom_no = 3; this_atom->atom_mass = 6.9675; strcpy(this_atom->symbol, "Li");  break;
        case 4: this_atom->atom_no = 4; this_atom->atom_mass = 9.01218; strcpy(this_atom->symbol, "Be");  break;
        case 5: this_atom->atom_no = 5; this_atom->atom_mass = 10.8135; strcpy(this_atom->symbol, "B");  break;
        case 6: this_atom->atom_no = 6; this_atom->atom_mass = 12.0106; strcpy(this_atom->symbol, "C");  break;
        case 7: this_atom->atom_no = 7; this_atom->atom_mass = 14.0069; strcpy(this_atom->symbol, "N");  break;
        case 8: this_atom->atom_no = 8; this_atom->atom_mass = 15.9994; strcpy(this_atom->symbol, "O");  break;
        case 9: this_atom->atom_no = 9; this_atom->atom_mass = 18.9984; strcpy(this_atom->symbol, "F");  break;
        case 10: this_atom->atom_no = 10; this_atom->atom_mass = 20.17976; strcpy(this_atom->symbol, "Ne");  break;
        case 11: this_atom->atom_no = 11; this_atom->atom_mass = 22.98977; strcpy(this_atom->symbol, "Na");  break;
        case 12: this_atom->atom_no = 12; this_atom->atom_mass = 24.3055; strcpy(this_atom->symbol, "Mg");  break;
        case 13: this_atom->atom_no = 13; this_atom->atom_mass = 26.98154; strcpy(this_atom->symbol, "Al");  break;
        case 14: this_atom->atom_no = 14; this_atom->atom_mass = 28.085; strcpy(this_atom->symbol, "Si");  break;
        case 15: this_atom->atom_no = 15; this_atom->atom_mass = 30.97376; strcpy(this_atom->symbol, "P");  break;
        case 16: this_atom->atom_no = 16; this_atom->atom_mass = 32.0675; strcpy(this_atom->symbol, "S");  break;
        case 17: this_atom->atom_no = 17; this_atom->atom_mass = 35.4515; strcpy(this_atom->symbol, "Cl");  break;
        case 18: this_atom->atom_no = 18; this_atom->atom_mass = 39.9481; strcpy(this_atom->symbol, "Ar");  break;
        case 19: this_atom->atom_no = 19; this_atom->atom_mass = 39.09831; strcpy(this_atom->symbol, "K");  break;
        case 20: this_atom->atom_no = 20; this_atom->atom_mass = 40.0784; strcpy(this_atom->symbol, "Ca");  break;
        case 21: this_atom->atom_no = 21; this_atom->atom_mass = 44.95591; strcpy(this_atom->symbol, "Sc");  break;
        case 22: this_atom->atom_no = 22; this_atom->atom_mass = 47.8671; strcpy(this_atom->symbol, "Ti");  break;
        case 23: this_atom->atom_no = 23; this_atom->atom_mass = 50.94151; strcpy(this_atom->symbol, "V");  break;
        case 24: this_atom->atom_no = 24; this_atom->atom_mass = 51.99616; strcpy(this_atom->symbol, "Cr");  break;
        case 25: this_atom->atom_no = 25; this_atom->atom_mass = 54.93804; strcpy(this_atom->symbol, "Mn");  break;
        case 26: this_atom->atom_no = 26; this_atom->atom_mass = 55.8452; strcpy(this_atom->symbol, "Fe");  break;
        case 27: this_atom->atom_no = 27; this_atom->atom_mass = 58.93319; strcpy(this_atom->symbol, "Co");  break;
        case 28: this_atom->atom_no = 28; this_atom->atom_mass = 58.69344; strcpy(this_atom->symbol, "Ni");  break;
        case 29: this_atom->atom_no = 29; this_atom->atom_mass = 63.5463; strcpy(this_atom->symbol, "Cu");  break;
        case 30: this_atom->atom_no = 30; this_atom->atom_mass = 65.382; strcpy(this_atom->symbol, "Zn");  break;
        case 31: this_atom->atom_no = 31; this_atom->atom_mass = 69.7231; strcpy(this_atom->symbol, "Ga");  break;
        case 32: this_atom->atom_no = 32; this_atom->atom_mass = 72.6308; strcpy(this_atom->symbol, "Ge");  break;
        case 33: this_atom->atom_no = 33; this_atom->atom_mass = 74.9216; strcpy(this_atom->symbol, "As");  break;
        case 34: this_atom->atom_no = 34; this_atom->atom_mass = 78.9718; strcpy(this_atom->symbol, "Se");  break;
        case 35: this_atom->atom_no = 35; this_atom->atom_mass = 79.904; strcpy(this_atom->symbol, "Br");  break;
        case 36: this_atom->atom_no = 36; this_atom->atom_mass = 83.7982; strcpy(this_atom->symbol, "Kr");  break;
        case 37: this_atom->atom_no = 37; this_atom->atom_mass = 85.46783; strcpy(this_atom->symbol, "Rb");  break;
        case 38: this_atom->atom_no = 38; this_atom->atom_mass = 87.621; strcpy(this_atom->symbol, "Sr");  break;
        case 39: this_atom->atom_no = 39; this_atom->atom_mass = 88.90584; strcpy(this_atom->symbol, "Y");  break;
        case 40: this_atom->atom_no = 40; this_atom->atom_mass = 91.2242; strcpy(this_atom->symbol, "Zr");  break;
        case 41: this_atom->atom_no = 41; this_atom->atom_mass = 92.90637; strcpy(this_atom->symbol, "Nb");  break;
        case 42: this_atom->atom_no = 42; this_atom->atom_mass = 95.951; strcpy(this_atom->symbol, "Mo");  break;
        case 43: this_atom->atom_no = 43; this_atom->atom_mass = 98.0; strcpy(this_atom->symbol, "Tc");  break;
        case 44: this_atom->atom_no = 44; this_atom->atom_mass = 101.072; strcpy(this_atom->symbol, "Ru");  break;
        case 45: this_atom->atom_no = 45; this_atom->atom_mass = 102.9055; strcpy(this_atom->symbol, "Rh");  break;
        case 46: this_atom->atom_no = 46; this_atom->atom_mass = 106.421; strcpy(this_atom->symbol, "Pd");  break;
        case 47: this_atom->atom_no = 47; this_atom->atom_mass = 107.86822; strcpy(this_atom->symbol, "Ag");  break;
        case 48: this_atom->atom_no = 48; this_atom->atom_mass = 112.4144; strcpy(this_atom->symbol, "Cd");  break;
        case 49: this_atom->atom_no = 49; this_atom->atom_mass = 114.8181; strcpy(this_atom->symbol, "In");  break;
        case 50: this_atom->atom_no = 50; this_atom->atom_mass = 118.7107; strcpy(this_atom->symbol, "Sn");  break;
        case 51: this_atom->atom_no = 51; this_atom->atom_mass = 121.7601; strcpy(this_atom->symbol, "Sb");  break;
        case 52: this_atom->atom_no = 52; this_atom->atom_mass = 127.603; strcpy(this_atom->symbol, "Te");  break;
        case 53: this_atom->atom_no = 53; this_atom->atom_mass = 126.90447; strcpy(this_atom->symbol, "I");  break;
        case 54: this_atom->atom_no = 54; this_atom->atom_mass = 131.2936; strcpy(this_atom->symbol, "Xe");  break;
        case 55: this_atom->atom_no = 55; this_atom->atom_mass = 132.90545; strcpy(this_atom->symbol, "Cs");  break;
        case 56: this_atom->atom_no = 56; this_atom->atom_mass = 137.3277; strcpy(this_atom->symbol, "Ba");  break;
        case 57: this_atom->atom_no = 57; this_atom->atom_mass = 138.90548; strcpy(this_atom->symbol, "La");  break;
        case 58: this_atom->atom_no = 58; this_atom->atom_mass = 140.1161; strcpy(this_atom->symbol, "Ce");  break;
        case 59: this_atom->atom_no = 59; this_atom->atom_mass = 140.90766; strcpy(this_atom->symbol, "Pr");  break;
        case 60: this_atom->atom_no = 60; this_atom->atom_mass = 144.2423; strcpy(this_atom->symbol, "Nd");  break;
        case 61: this_atom->atom_no = 61; this_atom->atom_mass = 145.0; strcpy(this_atom->symbol, "Pm");  break;
        case 62: this_atom->atom_no = 62; this_atom->atom_mass = 150.362; strcpy(this_atom->symbol, "Sm");  break;
        case 63: this_atom->atom_no = 63; this_atom->atom_mass = 151.9641; strcpy(this_atom->symbol, "Eu");  break;
        case 64: this_atom->atom_no = 64; this_atom->atom_mass = 157.253; strcpy(this_atom->symbol, "Gd");  break;
        case 65: this_atom->atom_no = 65; this_atom->atom_mass = 158.92535; strcpy(this_atom->symbol, "Tb");  break;
        case 66: this_atom->atom_no = 66; this_atom->atom_mass = 162.5001; strcpy(this_atom->symbol, "Dy");  break;
        case 67: this_atom->atom_no = 67; this_atom->atom_mass = 164.93033; strcpy(this_atom->symbol, "Ho");  break;
        case 68: this_atom->atom_no = 68; this_atom->atom_mass = 167.2593; strcpy(this_atom->symbol, "Er");  break;
        case 69: this_atom->atom_no = 69; this_atom->atom_mass = 168.93422; strcpy(this_atom->symbol, "Tm");  break;
        case 70: this_atom->atom_no = 70; this_atom->atom_mass = 173.0545; strcpy(this_atom->symbol, "Yb");  break;
        case 71: this_atom->atom_no = 71; this_atom->atom_mass = 174.96681; strcpy(this_atom->symbol, "Lu");  break;
        case 72: this_atom->atom_no = 72; this_atom->atom_mass = 178.492; strcpy(this_atom->symbol, "Hf");  break;
        case 73: this_atom->atom_no = 73; this_atom->atom_mass = 180.94788; strcpy(this_atom->symbol, "Ta");  break;
        case 74: this_atom->atom_no = 74; this_atom->atom_mass = 183.841; strcpy(this_atom->symbol, "W");  break;
        case 75: this_atom->atom_no = 75; this_atom->atom_mass = 186.2071; strcpy(this_atom->symbol, "Re");  break;
        case 76: this_atom->atom_no = 76; this_atom->atom_mass = 190.233; strcpy(this_atom->symbol, "Os");  break;
        case 77: this_atom->atom_no = 77; this_atom->atom_mass = 192.2173; strcpy(this_atom->symbol, "Ir");  break;
        case 78: this_atom->atom_no = 78; this_atom->atom_mass = 195.0849; strcpy(this_atom->symbol, "Pt");  break;
        case 79: this_atom->atom_no = 79; this_atom->atom_mass = 196.96657; strcpy(this_atom->symbol, "Au");  break;
        case 80: this_atom->atom_no = 80; this_atom->atom_mass = 200.5923; strcpy(this_atom->symbol, "Hg");  break;
        case 81: this_atom->atom_no = 81; this_atom->atom_mass = 204.384; strcpy(this_atom->symbol, "Tl");  break;
        case 82: this_atom->atom_no = 82; this_atom->atom_mass = 207.21; strcpy(this_atom->symbol, "Pb");  break;
        case 83: this_atom->atom_no = 83; this_atom->atom_mass = 208.9804; strcpy(this_atom->symbol, "Bi");  break;
        case 84: this_atom->atom_no = 84; this_atom->atom_mass = 209.0; strcpy(this_atom->symbol, "Po");  break;
        case 85: this_atom->atom_no = 85; this_atom->atom_mass = 210.0; strcpy(this_atom->symbol, "At");  break;
        case 86: this_atom->atom_no = 86; this_atom->atom_mass = 222.0; strcpy(this_atom->symbol, "Rn");  break;
        case 87: this_atom->atom_no = 87; this_atom->atom_mass = 223.0; strcpy(this_atom->symbol, "Fr");  break;
        case 88: this_atom->atom_no = 88; this_atom->atom_mass = 226.0; strcpy(this_atom->symbol, "Ra");  break;
        case 89: this_atom->atom_no = 89; this_atom->atom_mass = 227.0; strcpy(this_atom->symbol, "Ac");  break;
        case 90: this_atom->atom_no = 90; this_atom->atom_mass = 232.03774; strcpy(this_atom->symbol, "Th");  break;
        case 91: this_atom->atom_no = 91; this_atom->atom_mass = 231.03588; strcpy(this_atom->symbol, "Pa");  break;
        case 92: this_atom->atom_no = 92; this_atom->atom_mass = 238.02891; strcpy(this_atom->symbol, "U");  break;
        case 93: this_atom->atom_no = 93; this_atom->atom_mass = 237.0; strcpy(this_atom->symbol, "Np");  break;
        case 94: this_atom->atom_no = 94; this_atom->atom_mass = 244.1; strcpy(this_atom->symbol, "Pu");  break;
        case 95: this_atom->atom_no = 95; this_atom->atom_mass = 243.1; strcpy(this_atom->symbol, "Am");  break;
        case 96: this_atom->atom_no = 96; this_atom->atom_mass = 247.1; strcpy(this_atom->symbol, "Cm");  break;
        case 97: this_atom->atom_no = 97; this_atom->atom_mass = 247.1; strcpy(this_atom->symbol, "Bk");  break;
        case 98: this_atom->atom_no = 98; this_atom->atom_mass = 251.1; strcpy(this_atom->symbol, "Cf");  break;
        case 99: this_atom->atom_no = 99; this_atom->atom_mass = 254.1; strcpy(this_atom->symbol, "Es");  break;
        case 100: this_atom->atom_no = 100; this_atom->atom_mass = 257.1; strcpy(this_atom->symbol, "Fm");  break;
        case 101: this_atom->atom_no = 101; this_atom->atom_mass = 258.0; strcpy(this_atom->symbol, "Md");  break;
        case 102: this_atom->atom_no = 102; this_atom->atom_mass = 259.0; strcpy(this_atom->symbol, "No");  break;
        case 103: this_atom->atom_no = 103; this_atom->atom_mass = 262.0; strcpy(this_atom->symbol, "Lr");  break;
        case 104: this_atom->atom_no = 104; this_atom->atom_mass = 267.0; strcpy(this_atom->symbol, "Rf");  break;
        case 105: this_atom->atom_no = 105; this_atom->atom_mass = 268.0; strcpy(this_atom->symbol, "Db");  break;
        case 106: this_atom->atom_no = 106; this_atom->atom_mass = 271.0; strcpy(this_atom->symbol, "Sg");  break;
        case 107: this_atom->atom_no = 107; this_atom->atom_mass = 272.0; strcpy(this_atom->symbol, "Bh");  break;
        case 108: this_atom->atom_no = 108; this_atom->atom_mass = 270.0; strcpy(this_atom->symbol, "Hs");  break;
        case 109: this_atom->atom_no = 109; this_atom->atom_mass = 276.0; strcpy(this_atom->symbol, "Mt");  break;
        case 110: this_atom->atom_no = 110; this_atom->atom_mass = 281.0; strcpy(this_atom->symbol, "Ds");  break;
        case 111: this_atom->atom_no = 111; this_atom->atom_mass = 280.0; strcpy(this_atom->symbol, "Rg");  break;
        case 112: this_atom->atom_no = 112; this_atom->atom_mass = 285.0; strcpy(this_atom->symbol, "Cn");  break;
        case 113: this_atom->atom_no = 113; this_atom->atom_mass = 284.0; strcpy(this_atom->symbol, "Nh");  break;
        case 114: this_atom->atom_no = 114; this_atom->atom_mass = 289.0; strcpy(this_atom->symbol, "Fl");  break;
        case 115: this_atom->atom_no = 115; this_atom->atom_mass = 288.0; strcpy(this_atom->symbol, "Mc");  break;
        case 116: this_atom->atom_no = 116; this_atom->atom_mass = 293.0; strcpy(this_atom->symbol, "Lv");  break;
        case 117: this_atom->atom_no = 117; this_atom->atom_mass = 292.0; strcpy(this_atom->symbol, "Ts");  break;
        case 118: this_atom->atom_no = 118; this_atom->atom_mass = 294.0; strcpy(this_atom->symbol, "Og");  break;
        default: return 1;
    }
    return 0;
}

int element_to_no(char *element)
{
    int atom_no = -1;
    atominfo this_atom;
    for (int test_no = 1; test_no < 119; test_no++)
    {
        no_to_atominfo(test_no, &this_atom);
        if (strcmp(element, this_atom.symbol) == 0)
        {
            atom_no = this_atom.atom_no;
        }
    }
    return atom_no;
}

int no_to_element(int atom_no, char* element)
{
    atominfo this_atom;
    if (no_to_atominfo(atom_no, &this_atom) != 0)
    {
        return -1;
    }
    else
    {
        strcpy(element, this_atom.symbol);
        return 0;
    }
}

float no_to_mass(int atom_no)
{
    atominfo this_atom;
    if (no_to_atominfo(atom_no, &this_atom) != 0)
    {
        return -1;
    }
    else
    {
        return this_atom.atom_mass;
    }
}

float pbc_dist(float coord_1[3], float coord_2[3], float pbc[3][3])
{
    float dist;
    if ((pbc[0][1] != 0) || (pbc[0][2] != 0) || (pbc[1][0] != 0) || (pbc[1][2] != 0) || (pbc[2][0] != 0) || (pbc[2][1] != 0))
    {
        dist = pbc_dist_triclinic(coord_1, coord_2, pbc);
    }
    else
    {
        dist = pbc_dist_orthogonal(coord_1, coord_2, pbc);
    }
    return dist;
}

float pbc_dist_triclinic(float coord_1[3], float coord_2[3], float pbc[3][3])
{
        printf("Non-orthogonal PBCs are not supported yet. Sorry!\n");
        exit(1);
}

float pbc_dist_orthogonal(float coord_1[3], float coord_2[3], float pbc[3][3])
{
    float diff[3];
    float dist = 0;
    for (int i = 0; i < 3; i++)
    {
        diff[i] = coord_2[i] - coord_1[i];
        if (diff[i] < 0) 
            diff[i] *= -1;
        diff[i] = diff[i] - floor(diff[i] / pbc[i][i]) * pbc[i][i];
        if (diff[i] > pbc[i][i]/2)
        {
            diff[i] = pbc[i][i] - diff[i];
        }

        dist += diff[i] * diff[i];
    }
    dist = sqrt(dist);
    return dist;
}

// Finds closest neighbor of type atom_neighbor in trajectory defined by traj_frame and atom to coordinates coord; returns -1 if nothing is found
int nextneighbor_in_traj(float coord[3], int atom_no, float traj_frame[atom_no][3], int atom[atom_no], int atom_neighbor, float pbc[3][3])
{
    int index_neighbor = -1;
    float min_distance = pbc[0][0] + pbc[1][1] + pbc[2][2]; // larger distance than any two atoms can have between them
    for (int i = 0; i < atom_no; i++)
    {
        if (atom[i] == atom_neighbor || atom_neighbor == 0)
        {
            float dist = pbc_dist(coord, traj_frame[i], pbc);
            if (dist < min_distance)
            {
                min_distance = dist;
                index_neighbor = i;
            }
        }
        else 
            continue;
    }

    return index_neighbor;
}
