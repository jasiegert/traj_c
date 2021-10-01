#include <string.h>
#include <stdio.h>

int element_to_no(char *element)
{
    int atom_no;
    if (strcmp(element, "H") == 0) atom_no = 1;
    else if (strcmp(element, "Li") == 0) atom_no = 3;
    else if (strcmp(element, "O") == 0) atom_no = 8;
    else if (strcmp(element, "C") == 0) atom_no = 6;
    else if (strcmp(element, "N") == 0) atom_no = 7;
    else if (strcmp(element, "O") == 0) atom_no = 8;
    else if (strcmp(element, "P") == 0) atom_no = 15;
    else if (strcmp(element, "S") == 0) atom_no = 16;
    else if (strcmp(element, "Cl") == 0) atom_no = 17;
    else if (strcmp(element, "Cs") == 0) atom_no = 55;
    else return -1;

    return atom_no;
}

int no_to_element(int atom_no, char* element)
{
    switch (atom_no)
    {
        case 1: strcpy(element, "H"); break;
        case 3: strcpy(element, "Li"); break;
        case 6: strcpy(element, "C"); break;
        case 7: strcpy(element, "N"); break;
        case 8: strcpy(element, "O"); break;
        case 15: strcpy(element, "P"); break;
        case 16: strcpy(element, "S"); break;
        case 17: strcpy(element, "Cl"); break;
        case 55: strcpy(element, "Cs"); break;
        default: return -1;
    }
    return 0;
}

float no_to_mass(int atom_no)
{
    float atom_mass;
    switch (atom_no)
    {
        case 1: atom_mass = 1.0; break;
        case 3: atom_mass = 0.0; break;
        case 6: atom_mass = 12.0; break;
        case 7: atom_mass = 14.0; break;
        case 8: atom_mass = 16.0; break;
        case 15: atom_mass = 31.0; break;
        case 16: atom_mass = 32.0; break;
        case 17: atom_mass = 35.45; break;
        case 55: atom_mass = 0.0; break;
        default: return -1.0;
    }
//    if (atom_mass == 0)
//    {
//        char atom_label[3];
//        no_to_element(atom_no, atom_label);
//        printf("Warning: atom mass of %s is %f.", atom_label, atom_mass);
 //   }
    return atom_mass;
}
