#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "umfpack_interface.h"
#include "iterative_method.h"
#include "write.h"
#include "coarse_grid.h"
#include "fine_grid.h"
#include "convert_coord.h"
#include "residu.h"
#include "printVector.h"
#include "matrixOperations.h"

int two_grid(double *u, double *res, int start, int it_pre, int it_post, int *ia, int *ja, double *a, double *b, double *born_m, int *born_coord, int n, int m, int nR, double L)
/*
    But
    ===
    Résoudre le système linéaire n x n
                        
                             Au = b

    en utilisant la 'two-grid method' qui consiste à alterner des méthodes
    de lissage (Gauss)Seidel symétrisée) avec des corrections 'coarse-grid'.
 
    Les étapes de résolutions sont les suivantes :
        (1) PRE-SMOOTHING :
            la méthode Gauss-Seidel -INFERIEURE- est d'abord appliquée
            pendant 'it_pre' itérations afin d'avoir une pemière
            estimation de la solution du système
            + amortit les hautes fréquences des erreurs
        (2) COARSE-GRID CORRECTION :
            utilise le solveur directe UMFPACK afin de résoudre le système
            sur la coarse-grid, càd. une grille avec un pas de
            discrétisation 2 fois plus grand : hR = 2*h
            + amortit les basses fréquences des erreurs
        (3) POST-SMOOTHING :
            la méthode Gauss-Seidel -SUPERIEURE- est d'abord appliquée
            pendant 'it_post' itérations afin d'avoir une pemière
            estimation de la solution du système
            + amortit les hautes fréquences des erreurs
            
    ALGORITHME : résoudre A1.u = b
        (a) appliquer 'it_pre' itérations de Gauss-Seidel inférieur  -> u1
        (b) restriction du résidu : r2 = R1.(b - A1.u1)              -> r2
        (c) résoudre au 'coarsest level' avec UMFPACK : A1.cR = r1   -> cR
        (d) prolonger la correction calculée : u1' = u1 + P.cR       -> u1'
        (c) appliquer 'it_post' itérations de Gauss-Seidel supérieur -> u
        

    Arguments
    =========
    u           (output) - pointeur vers le tbaleau 'u' du vecteur solution
    res         (output) - pointeur vers le tableau 'res' du vecteur résidu
                           conteant les résidu de chaque étape et itération
    start       (input) - indice à partir duquel on commence à stocker
                           les résidus dans le tableau 'res'
    it_pre      (input) - nombre d'itérations pour le pre-smoothing
    it_post     (input) - nombre d'itérations pour le post-smoothing
    ia          (input) - pointeur vers le tableau 'ia' de la matrice A
    ja          (input) - pointeur vers le tableau 'ja' de la matrice A
    a           (input) - pointeur vers le tableau 'a' de la matrice A
    b           (input) - pointeur vers le tableau 'b'
    born_m      (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires en mètres
    born_coord  (input) - tableau contenant les coordonées des domaines
                          des parties recatngulaires
    n           (output) - pointeur vers le nombre d'inconus dans le système
    m           (input)  - nombre de points par direction dans la grille
    nR          (input) - nombre de parties rectangulaires
    L           (input)  - dimensions LxL de la grille (en mètres)
    
    => renvoit 0 si fonctionnement normal et 1 sinon
*/
{
    /* ----- DECLARATION DES VARIABLES ----- */
    // nombre d'inconnues sur la coarse grid
    int n_R;
    // Matrice A réduite de la coarse grid au format CSR
    int *iaR, *jaR;
    double *aR, *bR;
    // discrétisation correspondant à la coarse grid
    int mR = (m+1)/2;
    // incrément
    int i;
    // coordonnées des zones rectangulaires sur la coarse grid
    int *born_coordR;
    // matrice du résidu
    double *resMatP; // fine grid
    double *resMatR; // coarse grid
    // vecteur correction : Ac = r
    double *corrR; // coarse grid
    double *corrP; // fine grid
    
    /* ----- INITIALISATION COARSE GRID ----- */
    born_coordR  = malloc((nR * 4) * sizeof(int));
    if (born_coordR == NULL) {
        printf("\n ERREUR : pas de mémoire pour le vecteur des zones rect sur la coarse grid\n\n");
           return 1;
    }
    // on convertit les bornes (en mètres) en coordonnées pour la coarse grid
    createCoord(nR, L, mR, born_m, born_coordR);
    // on génère le problème -> on détermine la matrice A au format CSR
    // et le nombre d'inconnues pour la COARSE GRID
    if (prob(mR, &n_R, L, born_coordR, &iaR, &jaR, &aR, &bR, nR)) { return 1; }
    
    /* ----- ALOCATION DES TABLEAUX ----- */
    resMatR = malloc(n_R * sizeof(double));
    resMatP = malloc(n * sizeof(double));
    corrR = malloc(n_R * sizeof(double));
    corrP = calloc(n, sizeof(double));
    if (resMatR == NULL || resMatP == NULL || corrR == NULL || corrP == NULL) {
        printf("\n ERREUR : pas de mémoire pour vecteur résidu et correction\n\n");
           return 1;
    }
    /* ----- TWO-GRID ----- */
    // PRE-SMOOTHING
    for(i = 0; i < it_pre; i++) {
        gaussSeidelSymetrique(ia, ja, a, b, n, u);
        res[start + i] = residu(ia, ja, a, b, u, n);
    }
    // calcul du vecteur résidu 'resMatP'
    if(residuMatrix(ia, ja, a, b, u, n, resMatP)) { return 1; }
    // restriction du vecteur résidu 'resMatP' => 'resMatR'
    restriction (m, n_R, born_coord, born_coordR, nR, resMatR, resMatP, L);
    printVector(resMatP, L, born_coord, m, n, nR, "resMatP.txt");
    printVector(resMatR, L, born_coordR, mR, n_R, nR, "resMatR.txt");
    // calcul du vecteur correction 'corrR' sur la coarsest grid
    // solveur directe UMFPack résout le problème: Ac = r
    double t1 = mytimer(); // temps de solution
    if( solve_umfpack(n_R, iaR, jaR, aR, resMatR, corrR) )
        return 1;
    double t2 = mytimer(); // temps de solution
    printf("\nTemps de solution (CPU): %5.1f sec\n",t2-t1);
    // prolongation du vecteur correction 'corrR' => 'corrP'
    prolongation(m, n, L, born_coord, corrR, corrP, nR);
    // calcul du vecteur u amélioré par corrP
    add(u, corrP, n, 0);
    // calcul du résidu final
    res[start+it_pre] = residu(ia, ja, a, b, u, n);
    // POST-SMOOTHING
    for(i = 0; i < it_post; i++) {
        gaussSeidelSymetrique(ia, ja, a, b, n, u);
        res[start + it_pre + 1 +  i] = residu(ia, ja, a, b, u, n);
    }

    /* ----- LIBERER LA MEMOIRE DES VARIABLES GENERALES ----- */
    free(iaR); free(jaR); free(aR); free(bR); free(born_coordR);
    free(resMatR); free(resMatP); free(corrP); free(corrR);
    
    return 0;
}
