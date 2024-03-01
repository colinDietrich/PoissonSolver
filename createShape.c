#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "convert_coord.h"

int createCroix(double **born_m, int **born_coord, int *nR, double L, int m)
/*
    BUT
    ===
    cette fonction permet de créer un ensemble de parties rectangulaires
    de sorte à ce que la membrane aie une forme de croix
 
    ARGUMENTS
    =========
    born_m      (output) - pointeur vers le tableau contenant les
                           positions (en mètre) des domaines des parties
                           recatngulaires
    born_coord  (output) - pointeur vers le tableau contenant les
                           coordonées des domaines des parties
                           recatngulaires
    nR      (input) - pointeur vers le nombre de parties rectangulaire
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                   du domaine de la membrane
 */
{
    *nR = 4;
    *born_m = malloc((*nR * 4) * sizeof(double));
    *born_coord  = malloc((*nR * 4) * sizeof(int));
    // on vérifie si l'allocation est réussite
    if (*born_coord == NULL || *born_m == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour convertir en coordonnées\n\n");
        return 1;
    }
    (*born_m)[0] = 0.0; (*born_m)[1] = 0.6;
    (*born_m)[2] = 0.0; (*born_m)[3] = 0.6;
    (*born_m)[4] = 1.4; (*born_m)[5] = 2.0;
    (*born_m)[6] = 0.0; (*born_m)[7] = 0.6;
    (*born_m)[8] = 0.0; (*born_m)[9] = 0.6;
    (*born_m)[10] = 1.4; (*born_m)[11] = 2.0;
    (*born_m)[12] = 1.4; (*born_m)[13] = 2.0;
    (*born_m)[14] = 1.4; (*born_m)[15] = 2.0;
    // on convertit les bornes en coordonnées
    createCoord(*nR, L, m, *born_m, *born_coord);
    return 0;
}


int createP8(double **born_m, int **born_coord, int *nR, double L, int m)
/*
    BUT
    ===
    cette fonction permet de créer un ensemble de parties rectangulaires
    correspondant au projet numéro 8
 
    ARGUMENTS
    =========
    born_m      (output) - pointeur vers le tableau contenant les
                           positions (en mètre) des domaines des parties
                           recatngulaires
    born_coord  (output) - pointeur vers le tableau contenant les
                           coordonées des domaines des parties
                           recatngulaires
    nR      (input) - pointeur vers le nombre de parties rectangulaire
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                   du domaine de la membrane
 */
{
    *nR = 1;
    *born_m = malloc((*nR * 4) * sizeof(double));
    *born_coord  = malloc((*nR * 4) * sizeof(int));
    // on vérifie si l'allocation est réussite
    if (*born_coord == NULL || *born_m == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour convertir en coordonnées\n\n");
        return 1;
    }
    (*born_m)[0] = 0.4; (*born_m)[1] = 1.6;
    (*born_m)[2] = 1.8; (*born_m)[3] = 2.0;
    // on convertit les bornes en coordonnées
    createCoord(*nR, L, m, *born_m, *born_coord);
    return 0;
}

int createCroixInverse(double **born_m, int **born_coord, int *nR, double L, int m)
/*
    BUT
    ===
    cette fonction permet de créer un ensemble de parties rectangulaires
    de sorte à ce que la membrane aie une forme de croix
 
    ARGUMENTS
    =========
    born_m      (output) - pointeur vers le tableau contenant les
                           positions (en mètre) des domaines des parties
                           recatngulaires
    born_coord  (output) - pointeur vers le tableau contenant les
                           coordonées des domaines des parties
                           recatngulaires
    nR      (input) - pointeur vers le nombre de parties rectangulaire
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                   du domaine de la membrane
 */
{
    *nR = 3;
    *born_m = malloc((*nR * 4) * sizeof(double));
    *born_coord  = malloc((*nR * 4) * sizeof(int));
    // on vérifie si l'allocation est réussite
    if (*born_coord == NULL || *born_m == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour convertir en coordonnées\n\n");
        return 1;
    }
    (*born_m)[0] = 0.8; (*born_m)[1] = 1.2;
    (*born_m)[2] = 0.0; (*born_m)[3] = 2.0;
    (*born_m)[4] = 0.0; (*born_m)[5] = 0.8;
    (*born_m)[6] = 0.8; (*born_m)[7] = 1.2;
    (*born_m)[8] = 1.2; (*born_m)[9] = 2.0;
    (*born_m)[10] = 0.8; (*born_m)[11] = 1.2;
    // on convertit les bornes en coordonnées
    createCoord(*nR, L, m, *born_m, *born_coord);
    return 0;
}

int createCircle(double **born_m, int **born_coord, int *nR, double L, int m)
/*
    BUT
    ===
    cette fonction permet de créer un ensemble de parties rectangulaires
    de sorte à ce que la membrane aie une forme de cercle
 
    ARGUMENTS
    =========
    born_m      (output) - pointeur vers le tableau contenant les
                           positions (en mètre) des domaines des parties
                           recatngulaires
    born_coord  (output) - pointeur vers le tableau contenant les
                           coordonées des domaines des parties
                           recatngulaires
    nR      (input) - pointeur vers le nombre de parties rectangulaire
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                   du domaine de la membrane
 */
{
    *nR = 48;
    *born_m = malloc((*nR * 4) * sizeof(double));
    *born_coord  = malloc((*nR * 4) * sizeof(int));
    // on vérifie si l'allocation est réussite
    if (*born_coord == NULL || *born_m == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour convertir en coordonnées\n\n");
        return 1;
    }
    
    // coté en bas à droite
    (*born_m)[0] = 1.2; (*born_m)[1] = 2.0; (*born_m)[2] = 0.0; (*born_m)[3] = 0.05;
    (*born_m)[4] = 1.35; (*born_m)[5] = 2.0; (*born_m)[6] = 0.05; (*born_m)[7] = 0.1;
    (*born_m)[8] = 1.45; (*born_m)[9] = 2.0; (*born_m)[10] = 0.1; (*born_m)[11] = 0.15;
    (*born_m)[12] = 1.55; (*born_m)[13] = 2.0; (*born_m)[14] = 0.15; (*born_m)[15] = 0.2;
    (*born_m)[16] = 1.6; (*born_m)[17] = 2.0; (*born_m)[18] = 0.2; (*born_m)[19] = 0.25;
    (*born_m)[20] = 1.65; (*born_m)[21] = 2.0; (*born_m)[22] = 0.25; (*born_m)[23] = 0.3;
    (*born_m)[24] = 1.7; (*born_m)[25] = 2.0; (*born_m)[26] = 0.3; (*born_m)[27] = 0.35;
    (*born_m)[28] = 1.75; (*born_m)[29] = 2.0; (*born_m)[30] = 0.35; (*born_m)[31] = 0.4;
    (*born_m)[32] = 1.8; (*born_m)[33] = 2.0; (*born_m)[34] = 0.4; (*born_m)[35] = 0.45;
    (*born_m)[36] = 1.85; (*born_m)[37] = 2.0; (*born_m)[38] = 0.45; (*born_m)[39] = 0.55;
    (*born_m)[40] = 1.9; (*born_m)[41] = 2.0; (*born_m)[42] = 0.55; (*born_m)[43] = 0.65;
    (*born_m)[44] = 1.95; (*born_m)[45] = 2.0; (*born_m)[46] = 0.65; (*born_m)[47] = 0.8;
    
    // coté en haut à droite
    (*born_m)[48] = 1.2; (*born_m)[49] = 2.0; (*born_m)[50] = 1.95; (*born_m)[51] = 2.0;
    (*born_m)[52] = 1.35; (*born_m)[53] = 2.0; (*born_m)[54] = 1.9; (*born_m)[55] = 1.95;
    (*born_m)[56] = 1.45; (*born_m)[57] = 2.0; (*born_m)[58] = 1.85; (*born_m)[59] = 1.9;
    (*born_m)[60] = 1.55; (*born_m)[61] = 2.0; (*born_m)[62] = 1.8; (*born_m)[63] = 1.85;
    (*born_m)[64] = 1.6; (*born_m)[65] = 2.0; (*born_m)[66] = 1.75; (*born_m)[67] = 1.8;
    (*born_m)[68] = 1.65; (*born_m)[69] = 2.0; (*born_m)[70] = 1.7; (*born_m)[71] = 1.75;
    (*born_m)[72] = 1.7; (*born_m)[73] = 2.0; (*born_m)[74] = 1.65; (*born_m)[75] = 1.7;
    (*born_m)[76] = 1.75; (*born_m)[77] = 2.0; (*born_m)[78] = 1.6; (*born_m)[79] = 1.65;
    (*born_m)[80] = 1.8; (*born_m)[81] = 2.0; (*born_m)[82] = 1.55; (*born_m)[83] = 1.6;
    (*born_m)[84] = 1.85; (*born_m)[85] = 2.0; (*born_m)[86] = 1.45; (*born_m)[87] = 1.55;
    (*born_m)[88] = 1.9; (*born_m)[89] = 2.0; (*born_m)[90] = 1.35; (*born_m)[91] = 1.45;
    (*born_m)[92] = 1.95; (*born_m)[93] = 2.0; (*born_m)[94] = 1.2; (*born_m)[95] = 1.35;
    
    // coté en bas à gauche
    (*born_m)[96] = 0.0; (*born_m)[97] = 0.05; (*born_m)[98] = 0.65; (*born_m)[99] = 0.8;
    (*born_m)[100] = 0.0; (*born_m)[101] = 0.1; (*born_m)[102] = 0.55; (*born_m)[103] = 0.65;
    (*born_m)[104] = 0.0; (*born_m)[105] = 0.15; (*born_m)[106] = 0.45; (*born_m)[107] = 0.55;
    (*born_m)[108] = 0.0; (*born_m)[109] = 0.2; (*born_m)[110] = 0.4; (*born_m)[111] = 0.45;
    (*born_m)[112] = 0.0; (*born_m)[113] = 0.25; (*born_m)[114] = 0.35; (*born_m)[115] = 0.4;
    (*born_m)[116] = 0.0; (*born_m)[117] = 0.3; (*born_m)[118] = 0.3; (*born_m)[119] = 0.35;
    (*born_m)[120] = 0.0; (*born_m)[121] = 0.35; (*born_m)[122] = 0.25; (*born_m)[123] = 0.3;
    (*born_m)[124] = 0.0; (*born_m)[125] = 0.4; (*born_m)[126] = 0.2; (*born_m)[127] = 0.25;
    (*born_m)[128] = 0.0; (*born_m)[129] = 0.45; (*born_m)[130] = 0.15; (*born_m)[131] = 0.2;
    (*born_m)[132] = 0.0; (*born_m)[133] = 0.55; (*born_m)[134] = 0.1; (*born_m)[135] = 0.15;
    (*born_m)[136] = 0.0; (*born_m)[137] = 0.65; (*born_m)[138] = 0.05; (*born_m)[139] = 0.1;
    (*born_m)[140] = 0.0; (*born_m)[141] = 0.8; (*born_m)[142] = 0.0; (*born_m)[143] = 0.05;
    
    // coté en haut à gauche
    (*born_m)[144] = 0.0; (*born_m)[145] = 0.05; (*born_m)[146] = 1.2; (*born_m)[147] = 1.35;
    (*born_m)[148] = 0.0; (*born_m)[149] = 0.1; (*born_m)[150] = 1.35; (*born_m)[151] = 1.45;
    (*born_m)[152] = 0.0; (*born_m)[153] = 0.15; (*born_m)[154] = 1.45; (*born_m)[155] = 1.55;
    (*born_m)[156] = 0.0; (*born_m)[157] = 0.2; (*born_m)[158] = 1.55; (*born_m)[159] = 1.6;
    (*born_m)[160] = 0.0; (*born_m)[161] = 0.25; (*born_m)[162] = 1.6; (*born_m)[163] = 1.65;
    (*born_m)[164] = 0.0; (*born_m)[165] = 0.3; (*born_m)[166] = 1.65; (*born_m)[167] = 1.7;
    (*born_m)[168] = 0.0; (*born_m)[169] = 0.35; (*born_m)[170] = 1.7; (*born_m)[171] = 1.75;
    (*born_m)[172] = 0.0; (*born_m)[173] = 0.4; (*born_m)[174] = 1.75; (*born_m)[175] = 1.8;
    (*born_m)[176] = 0.0; (*born_m)[177] = 0.45; (*born_m)[178] = 1.8; (*born_m)[179] = 1.85;
    (*born_m)[180] = 0.0; (*born_m)[181] = 0.55; (*born_m)[182] = 1.85; (*born_m)[183] = 1.9;
    (*born_m)[184] = 0.0; (*born_m)[185] = 0.65; (*born_m)[186] = 1.9; (*born_m)[187] = 1.95;
    (*born_m)[188] = 0.0; (*born_m)[189] = 0.8; (*born_m)[190] = 1.95; (*born_m)[191] = 2.0;
    
    // on convertit les bornes en coordonnées
    createCoord(*nR, L, m, *born_m, *born_coord);
    return 0;
}

int createP15(double **born_m, int **born_coord, int *nR, double L, int m)
/*
    BUT
    ===
    cette fonction permet de créer un ensemble de parties rectangulaires
    correspondant au projet numéro 8
 
    ARGUMENTS
    =========
    born_m      (output) - pointeur vers le tableau contenant les
                           positions (en mètre) des domaines des parties
                           recatngulaires
    born_coord  (output) - pointeur vers le tableau contenant les
                           coordonées des domaines des parties
                           recatngulaires
    nR      (input) - pointeur vers le nombre de parties rectangulaire
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                   du domaine de la membrane
 */
{
    *nR = 1;
    *born_m = malloc((*nR * 4) * sizeof(double));
    *born_coord  = malloc((*nR * 4) * sizeof(int));
    // on vérifie si l'allocation est réussite
    if (*born_coord == NULL || *born_m == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour convertir en coordonnées\n\n");
        return 1;
    }
    (*born_m)[0] = 1.0/3.0; (*born_m)[1] = 5.0/6.0;
    (*born_m)[2] = 1.0/6.0; (*born_m)[3] = 7.0/12.0;
    // on convertit les bornes en coordonnées
    createCoord(*nR, L, m, *born_m, *born_coord);
    return 0;
}
