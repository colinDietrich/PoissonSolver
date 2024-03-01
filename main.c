#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "time.h"
#include "check_discretisation.h"
#include "convert_coord.h"
#include "residu.h"
#include "createShape.h"
#include "write.h"
#include "printVector.h"
#include "two_grid.h"
#include "multiGrid.h"
#include "interface_primme.h"
#include "umfpack_interface.h"
#include "matrixOperations.h"

int main(int argc, char *argv[])
/*
    But
    ===
    main relie toutes les fonctions relatives au projet.
    Il existe 4 options :
     -> voir les descriptions de chacunes ci-dessous
    
    nb : la variable 'manuel' permet de configurer les parties rectangulaires
         manuellement si sa valeur est mise sur 1. Si elle vaut 0 des fonctions
         définies dans le fichier 'creatShape.c' permet directement de créer des
         formes prédéfinies.
         -> le mode manuel permet de mettre autant de parties rectangulaires
            que l'on souhaite. Pour cela, il faut préciser le nombre
            de parties rectangulaires avec la variable 'nR'
            + la variable 'born_m' permet de stocker les bornes [x1,x2] x [y1,y2]
              de chaque partie rectangulaire sous la forme :
              born_m[i] = x1, born_m[i+1] = x2, born_m[i+2] = y1, born_m[i+3] = y2
*/
{
    
    
    /* ----- DECLARATION DES VARIABLES ----- */

    // --- GENERALE ---
    // dimension de la membrane
    int m = 641;
    double L = 2;
    // tolérance admise pour définir le bord de la mebrane
    double tol = 1e-4;
    // nombre de parties rectangulaire
    int nR;
    // coordonnées des parties rectangulaires la membrane est fixée (u = 0)
    double *born_m;
    int *born_coord;
    // matrice de discrétisation A au format CSR
    int *ia, *ja;
    double *a, *b;
    // nombre d'inconnues dans le système A*u = (w^2)*u
    int n;
    
    // --- MULTI-GRID ---
    // nombre de niveaux pour le multi-grid
    int totalLevels = 2;
    // nombre d'itérations pour la multi-grid
    int iterations = 20;
    // nombre d'itérations pour le pre-smoothing
    int it_pre = 2;
    // nombre d'itérations pour le post-smoothing
    int it_post = 2;
    // facteur de relaxation
    double tau = 1.0;
    // type de cycle -> 1 = V-cycle et 2 = W-cycle
    int cycle = 1;
    // stockage des différentes valeurs de résidus
    int size = 1 + iterations;
    // incrément
    int i;
    // niveau actuel de l'algorithme multi-grid
    int level;
    // rafinement de la grille la plus grossière
    int mCoarse = (m + pow(2, totalLevels-1) - 1) / pow(2, totalLevels-1);
    // tableau des résidus de chaque itération
    double *res;
    // tableaux de pointeurs pour stocker les tableaux de tous les niveaux
    int **iaM, **jaM, **born_coordM, *nM;
    double **aM, **bM, **uM;
    // tableaux pour le préconditionneur multi-grid
    double **corrM, **resM, *sol;
    
    // --- OPTION D'EXECUTION ---
    // entree manuel ou non de la configuration de la membrane
    int manuel = 0;
    // numéro de l'option
    int option = 1;
    
    
    printf("\n****************************************");
    printf("\n*********   PROJET MATH-H401    ********");
    printf("\n*********    Dietrich Colin     ********");
    printf("\n********* matricule : 000474612 ********");
    printf("\n****************************************");
    
    /* -------------------------- */
    /* ----- INITIALISATION ----- */
    /* -------------------------- */
    
    // on configure la membrane
    if(!manuel) {
        if(createP8(&born_m, &born_coord, &nR, L, m)) { return 1; }
    }
    else {
        
        // --- EXEMPLE POUR PLUSIEURS RECTANGLES ---
        // nombre de parties rectangulaires
        nR = 1;
        // allouer la memoire pour les bornes
        born_m = malloc((nR * 4) * sizeof(double));
        born_coord  = malloc((nR * 4) * sizeof(int));
        // bornes en mètres
        born_m[0] = 0.2; born_m[1] = 0.3;
        born_m[2] = 0.0; born_m[3] = 0.7;
        // on convertit les bornes (en mètres) en coordonnées
        createCoord(nR, L, m, born_m, born_coord);
        // -----------------------------------------
        
    }
    
    // on vérifie que le problème est bien posé avec le paramètre m
    if(checkDiscretisation(m, L, born_m, nR, tol)) { return 1; }
    
    // on génère le problème -> on détermine la matrice A au format CSR
    // et le nombre d'inconnues
    if (prob(m, &n, L, born_coord, &ia, &ja, &a, &b, nR)) { return 1; }
    if(writeTabInt(n + 1, ia, "ia.txt")) { return 1; }
    if(writeTabInt(ia[n], ja, "ja.txt")) { return 1; }
    if(writeTabDouble(ia[n], a, "a.txt")) { return 1; }
    if(writeTabDouble(n, b, "b.txt")) { return 1; }
    
    printf("\n----------------------------------------");
    printf("\n--     INITIALISATION DU PROBLEME     --");
    printf("\n----------------------------------------\n");
    printf("\nrafinement : m = %d\n", m);
    printf("nombre d'inconnues : n = %d\n", n);
    printf("nombre d'éléments non nuls de la matrice A : nnz = %d\n\n", ia[n]);
    
    

    /* ----- OPTION 1 -----
     Cette option permet de visualiser l'impact de la correction coarse-grid
     par rapport aux itérations de lissage
     -> elle utilise la fonction 'two_grid' codée initialement pour réaliser
        la méthode two-grid
        /!\ Il est plutôt conseillé d'utiliser les option 2 ou 3 pour implémenter
            la méthode two-grid. Pour ce faire il suffit seulement de mettre la
            variable 'totalLevels' sur 2 et d'utiliser la fonction 'multiGrid'
            utilisée dans les options 1 et 2
     -> cette option a été utilisée pour réaliser la figure 8 du rapport
    ----------------------- */
    if(option == 1) {
        // stockage des différentes valeurs de résidus
        size = 1 + (it_pre + it_post + 1) * iterations;
        // incrément
        int i;
        
        /* ----- ALOCATION DES TABLEAUX ----- */
        sol = calloc(n, sizeof(double));
        res = malloc(size * sizeof(double));
        if (sol == NULL || res == NULL) {
            printf("\n ERREUR : pas de mémoire pour le vecteur solution et résidus\n\n");
              return 1;
        }
        
        /* ----- TWO-GRID METHOD ----- */
        // calcul du premier résidu
        res[0] = residu(ia, ja, a, b, sol, n);
        // itérations
        for(i = 0; i < iterations; i++) {
            int start = 1 + (it_pre + it_post) * i;
            if(two_grid(sol, res, start, it_pre, it_post, ia, ja, a, b, born_m, born_coord, n, m, nR, L)) { return 1; };
        }
        if(writeTabDouble(size, res, "res.txt")) { return 1; }
        free(sol); free(res);
    }
    
    /* ----- OPTION 2 -----
     - SOLVEUR MULTI-GRILLES -
     Cette option implémente la fonction 'multi-grid' pour résoudre l'équation de Poisson
     -> elle résout le système Au = b
     -> cette option ne permet pas d'implémenter le facteur de relaxation 'tau'
    ----------------------- */
    if(option == 2) {
        
        // on vérifie que le problème est bien posé avec les paramètres m et levels
        if(checkMultigrid(m, L, born_m, nR, tol, totalLevels)) { return 1; }
        
        /* ----- ALOCATION DES TABLEAUX ----- */
        res = malloc(size * sizeof(double));
        iaM = malloc(totalLevels * sizeof(int*));
        jaM = malloc(totalLevels * sizeof(int*));
        aM = malloc(totalLevels * sizeof(double*));
        bM = malloc(totalLevels * sizeof(double*));
        uM = malloc(totalLevels * sizeof(double*));
        born_coordM = malloc(totalLevels * sizeof(int*));
        nM = malloc(totalLevels * sizeof(int));
        if (iaM == NULL || jaM == NULL || aM == NULL || bM == NULL || uM == NULL || born_coordM == NULL || nM == NULL || res == NULL) {
            printf("\n ERREUR : pas de mémoire pour la méthode multi-grid\n\n");
              return 1;
        }
        /* ----- INITIALISATION ---- */
        if( init_multiGrid(totalLevels, nR, L, m, born_m, &iaM, &jaM, &aM, &bM, &uM, &born_coordM, &nM) ) { return 1; }
        /* ----- MULTI-GRID METHOD ----- */
        // calcul du preier résidu
        res[0] = residu(iaM[0], jaM[0], aM[0], bM[0], uM[0], nM[0]);
        // itérations
        for(i = 0; i < iterations; i++) {
            int j;
            for(j = 1; j < totalLevels; j++) {
                clearArray(uM[j], nM[j]);
            }
            level = -1; // réinitialisation du niveau
            if(multiGrid(it_pre, it_post, uM, iaM, jaM, aM, bM, born_coordM, nM, m, nR, L, mCoarse, level, cycle)) { return 1; }
            res[i + 1] = residu(iaM[0], jaM[0], aM[0], bM[0], uM[0], nM[0]);
        }
        if(writeTabDouble(size, res, "res.txt")) { return 1; }
        printVector(uM[0], L, born_coordM[0], m, nM[0], nR, "u.txt");
        
        /* ----- LIBERER MEMOIRE ----- */
        free(res); free(iaM); free(jaM); free(aM); free(bM); free(uM);
        free(born_coordM); free(nM);
    }
    
    /* ----- OPTION 3 -----
     - PRECONDITIONNEUR MULTI-GRILLES -
     Cette option implémente la fonction 'multi-grid' pour résoudre l'équation de Poisson
     -> elle résout le système Ax = b - Au et ajoute la correction sol += tau * x à chaque itération
     -> cette option permet d'implémenter le facteur de relaxation 'tau'
     -> cette option permet de calculer les valeurs propres min/max, le facteur de convergence
        asymptotique, le tau optimal et vérifie la stabilité de l'algorithme
    ----------------------- */
    if(option == 3) {
        /* ----- DECLARATION DES VARIABLES ----- */
        double tolMultigrid = 1e-20; // tolerance pour la convergence
        double res1, res2 = 2.0; // stockage des 2 derniers résidus
        double rho; // facteur de convergence
        double normU; // norme de u
        double h = L / (m-1); // pas de discrétisation h
        double l_max = 8.0/(h*h); // valeur propre max de A
        double stable; // vérification de la stabilité directe la multi-grid method
        double o_mach = 1.1e-16; // précision machine
        double l_min; // valeur propre min de A
        double tau_opt; // facteur de relaxation optimal
        int i = 0; // incrément
        
        // on vérifie que le problème est bien posé avec les paramètres m et levels
        if(checkMultigrid(m, L, born_m, nR, tol, totalLevels)) { return 1; }
        
        /* ----- ALOCATION DES TABLEAUX ----- */
        iaM = malloc(totalLevels * sizeof(int*));
        jaM = malloc(totalLevels * sizeof(int*));
        aM = malloc(totalLevels * sizeof(double*));
        resM = malloc(totalLevels * sizeof(double*));
        corrM = malloc(totalLevels * sizeof(double*));
        born_coordM = malloc(totalLevels * sizeof(int*));
        nM = malloc(totalLevels * sizeof(int));
        if (iaM == NULL || jaM == NULL || aM == NULL || resM == NULL || corrM == NULL || born_coordM == NULL || nM == NULL) {
            printf("\n ERREUR : pas de mémoire pour la méthode multi-grid\n\n");
              return 1;
        }
        /* ----- INITIALISATION ---- */
        if( init_multiGrid(totalLevels, nR, L, m, born_m, &iaM, &jaM, &aM, &resM, &corrM, &born_coordM, &nM) ) { return 1; }
        sol = calloc(nM[0], sizeof(double));
        if (sol == NULL) {
            printf("\n ERREUR : pas de mémoire pour vecteur sol\n\n");
              return 1;
        }
        /* ----- MULTI-GRID METHOD ----- */
        // calcul du preier résidu
        res1 = residu(iaM[0], jaM[0], aM[0], b, sol, nM[0]);
        appendText(res1, "resIter.txt");
        if(residuMatrix(iaM[0], jaM[0], aM[0], b, sol, nM[0], resM[0])) { return 1; }
        // itérations jusqu'à la convergence du résidu
        while(res2-res1 >= tolMultigrid) {
            level = -1; // réinitialisation du niveau
            int j;
            for(j = 0; j < totalLevels; j++) {
                clearArray(corrM[j], nM[j]);
            }
            if(multiGrid(it_pre, it_post, corrM, iaM, jaM, aM, resM, born_coordM, nM, m, nR, L, mCoarse, level, cycle)) { return 1; }
            for(j = 0; j < nM[0]; j++) {
                sol[j] += tau * corrM[0][j];
            }
            if(residuMatrix(iaM[0], jaM[0], aM[0], b, sol, nM[0], resM[0])) { return 1; }
            res2 = res1; // résidu de l'itération précédente
            res1 = residu(iaM[0], jaM[0], aM[0], b, sol, nM[0]); // nouveau résidu
            appendText(res1, "resIter.txt"); // enregistre les résidus successifs
            // on sélectionnne l'endroit où on veut calculer rho avec i
            if(i == 6) {
                rho = res1/res2;
            }
            i++;
        }
        printVector(sol, L, born_coordM[0], m, nM[0], nR, "u.txt");
        
        /* --- Calcul des paramètres de la question 4 --- */
        normU = norm(sol, n); // norme de sol
        stable = res1 / (l_max * normU);
        l_min = 1.0 - rho; // valeur propre max de A
        tau_opt = 2.0 / (1.0 + l_min); // facteur de relaxation optimal
        printf("--- VERIFICATION DE LA STABILITE ---\n");
        printf("    -> précision machine :   %e\n", o_mach);
        printf("    -> |r| / (|A|.|u|) :   %e\n", stable);
        printf("    -> facteur de convergence :   (%e/%e) = %e\n",res1, res2, rho);
        printf("    -> valeur propre minimale :   %e\n", l_min);
        printf("    -> facteur de relaxation optimal :   %e\n", tau_opt);
        printf("    -> l_max corrigée : %e\n", (1+rho)/tau);
        
        /* ----- LIBERER MEMOIRE ----- */
        free(iaM); free(jaM); free(aM); free(bM); free(corrM);
        free(born_coordM); free(nM); free(sol);
    }

    
    /* ----- OPTION 4 -----
     - PRECONDITIONNEUR MULTI-GRILLES -
     Cette option permet de résoudre le problème aux valeurs propre Au = w^2u
     avec le solveur PRIMME
     -> elle permet d'implémenter le préconditionneur multi-grid pour accéler
        la convergence de PRIMME
    ----------------------- */
    if(option == 4) {
        /* ----- DECLARATION DES VARIABLES ----- */
        // tableau des valeurs et vecteurs propres
        double *evalsPrimme, *evecsPrimme;
        // nombre de valeurs propres à calculer
        int nev = 1;
        // plus petite (0) / grande (1) valeur propre
        int largest = 0;
        // temps de calcul du solveur primme
        double tc1Primme, tc2Primme;
        
        /* ----- ALOCATION DES TABLEAUX ----- */
        evalsPrimme = malloc(nev * sizeof(double));
        evecsPrimme = malloc(nev * n * sizeof(double));
        if (evalsPrimme == NULL || evecsPrimme == NULL) {
            printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
            return 1;
        }
        
        /* ----- RESOLUTION AVEC PRIMME ----- */
        printf("\n----------------------------------------");
        printf("\n--           SOLVEUR PRIMME           --");
        printf("\n----------------------------------------\n");
        // primme - résolution -> recherche de valeur propre + vexteur propres correspondants
        tc1Primme = mytimer();
        if(primme(n, ia, ja, a, nev, evalsPrimme, evecsPrimme, largest, nR, L, m, born_m)) { return 1; }
        tc2Primme = mytimer();
        
        // affichage
        printf("\nTemps de solution (CPU): %e sec",tc2Primme-tc1Primme);
        printf("\nValeur propre numéro %d: %e\n",nev, evalsPrimme[0]);
        printVector(evecsPrimme, L, born_coord, m, n, nR, "evecs.txt");
        
    }
     
    /* ----- LIBERER LA MEMOIRE DES VARIABLES GENERALES ----- */
    free(ia); free(ja); free(a); free(b); free(born_m);
    free(born_coord);
    
    
    return 0;
}

