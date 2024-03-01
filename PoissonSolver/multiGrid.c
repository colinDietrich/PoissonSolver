#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "umfpack_interface.h"
#include "interface_primme.h"
#include "iterative_method.h"
#include "write.h"
#include "coarse_grid.h"
#include "fine_grid.h"
#include "convert_coord.h"
#include "residu.h"
#include "matrixOperations.h"
#include "printVector.h"

int multiGrid(int it_pre, int it_post, double **uM, int **iaM, int **jaM, double **aM, double **bM, int **born_coordM, int *nM, int m, int nR, double L, double mCoarse, int niv, int iter)
/*
    BUT
    ===
    cette fonction réalise la méthode multigrid
 
    ARGUMENTS
    =========
    it_pre          (input) - nombre d'itération de pré-smoothing
    it_post         (input) - nombre d'itération de post-smoothing
    uM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur u
                                des différents niveaux
    iaM             (output) - pointeur vers le tableau de pointeurs conteant les vecteur ia
                                des différents niveaux
    jaM             (output) - pointeur vers le tableau de pointeurs conteant les vecteur ja
                                des différents niveaux
    aM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur a
                                des différents niveaux
    bM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur b
                                des différents niveaux
    born_coordM     (output) - pointeur vers le tableau de pointeurs conteant les vecteur born_coord
                                des différents niveaux
    nM              (output) - pointeur vers le tableau de pointeurs conteant les nombres d'onconnues
                                des différents niveaux
    m               (input) - rafinement de la grille discrétisée
    nR              (input) - nombre de trous dans la membrane
    L               (input) - dimension en mètre de la membrane
    mCoarse         (input) - rafinement de la grille le plus grossière
    niv             (input) - niveau actuel de l'algorithme Multi-grid
    iter            (input) - si iter = 1  =>  cycle en V
                              si iter = 2  =>  cycle en W
 
    => renvoit 0 si fonctionne normalement et 0 sinon
 
    nb : pour observer le résidu à un endroit dans le code, il faut copier coller ce morceau :
         // visualisation du résidu
         if(niv == 0) {
             double *resBeforePre = malloc(nM[niv] * sizeof(double));
             if(residuMatrix(iaM[niv], jaM[niv], aM[niv], bM[niv], uM[niv], nM[niv], resBeforePre)) { return 1; }
             printVector(resBeforePre, L, born_coordM[niv], m, nM[niv], nR, "resBeforePre.txt");
             free(resBeforePre);
        }
*/
{
    
    /* ----- DECLARATION DES VARIABLES ----- */
    // rafinement de niveaux suivant
    double mR = (m+1)/2;
    // niveau local actuel
    niv ++;
    // vecteur résidu
    double *resAfterPre;
    // incrément
    int i;
    
    /* ----- MULTI-GRID ----- */    
    // PRE-SMOOTHING
    for(i = 0; i < it_pre; i++) {
        gaussSeidelSymetrique(iaM[niv], jaM[niv], aM[niv], bM[niv], nM[niv], uM[niv]);
    }
    
    // RESTRICTION du résidu
    resAfterPre = malloc(nM[niv] * sizeof(double));
    if(resAfterPre == NULL) {
        printf("\n ERREUR : pas de mémoire pour vecteur resAfterPre\n\n");
           return 1;
    }
    if(residuMatrix(iaM[niv], jaM[niv], aM[niv], bM[niv], uM[niv], nM[niv], resAfterPre)) { return 1; }
    // on génère la restriction du vecteur résidu 'resAfterPre'
    restriction (m, nM[niv+1], born_coordM[niv], born_coordM[niv+1], nR, bM[niv+1], resAfterPre, L);
        
    if(mR == mCoarse) {
        // SOLVEUR UMFPACK -> grille la plus grossière
        double t1 = mytimer();
        if( solve_umfpack(nM[niv+1], iaM[niv+1], jaM[niv+1], aM[niv+1], bM[niv+1], uM[niv+1]) )
            return 1;
        double t2 = mytimer();
        printf("\nTemps de solution (CPU): %5.1f sec\n",t2-t1);

    } else {
        for(i = 0; i < iter; i++) {
            multiGrid(it_pre, it_post, uM, iaM, jaM, aM, bM, born_coordM, nM, mR, nR, L, mCoarse, niv, iter);
        }
    }
    
    // PROLONGATION
    prolongation(m, nM[niv], L, born_coordM[niv], uM[niv+1], uM[niv], nR);
    
    // POST-SMOOTHING
    for(i = 0; i < it_post; i++) {
        gaussSeidelSymetrique(iaM[niv], jaM[niv], aM[niv], bM[niv], nM[niv], uM[niv]);
    }
    
    free(resAfterPre);
    return 0;
}

int init_multiGrid(int levels, int nR, int L, int m, double *born_m, int ***iaM, int ***jaM, double ***aM, double ***bM, double ***uM, int ***born_coordM, int **nM)
/*
    BUT
    ===
    cette fonction calcul et garde en mémoire toutes les matrices au format CSR des différents niveaux de
    l'lgorithme Multi-Grid ainsi que les coordonnées des trous correspondants.
    Il initialise également les vecteur b et u.
 
    ARGUMENTS
    =========
    levels          (input) - nombre de niveaux total de l'algorithme Multi-grid
    nR              (input) - nombre de trous dans la membrane
    L               (input) - dimension en mètre de la membrane
    m               (input) - rafinement de la grille discrétisée
    born_m          (input) - coordonnées en mètre des trous de la membrane
    iaM             (output) - pointeur vers le tableau de pointeurs conteant les vecteur ia
                            des différents niveaux
    jaM             (output) - pointeur vers le tableau de pointeurs conteant les vecteur ja
                            des différents niveaux
    aM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur a
                            des différents niveaux
    bM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur b
                            des différents niveaux
    uM              (output) - pointeur vers le tableau de pointeurs conteant les vecteur u
                            des différents niveaux
    born_coordM     (output) - pointeur vers le tableau de pointeurs conteant les vecteur born_coord
                            des différents niveaux
    nM              (output) - pointeur vers le tableau de pointeurs conteant les nombres d'onconnues
                            des différents niveaux
 
    => renvoit 0 si fonctionne normalement et 0 sinon
*/
{
    // incrément
    int i;
    for(i = 0; i < levels; i++) {
        // initialisation - born_coordM -
        (*born_coordM)[i] = malloc((nR * 4) * sizeof(int));
        if ((*born_coordM)[i] == NULL) {
            printf("\n ERREUR : pas de mémoire pour le vecteur des zones rect sur la fine grid\n\n");
               return 1;
        }
        createCoord(nR, L, m, born_m, (*born_coordM)[i]);
        // initialisation - iaM, jaM, aM, bM, nM -
        if (prob(m, &(*nM)[i], L, (*born_coordM)[i], &(*iaM)[i], &(*jaM)[i], &(*aM)[i], &(*bM)[i], nR)) { return 1; }
        // initialisation - uM -
        (*uM)[i] = calloc((*nM)[i], sizeof(double));
        if ((*uM)[i] == NULL) {
            printf("\n ERREUR : pas de mémoire pour vecteurs solutions\n\n");
               return 1;
        }
        m = (m+1)/2;
    }
    return 0;
}
