#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "unknown_calculator.h"
#include "conditionBord.h"
#include "dirichlet.h"

int prob(int m, int *n, double L, int *born_coord, int **ia, int **ja, double **a, double **b, int nR)
/*
    But
    ===
    Génerer le système linéaire n x n
                        
                             Au = b

    qui correspond à la disrétisation sur une grille cartesienne
    regulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,L] x [0,L]
           dx   dx       dy   dy

    avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
        u = exp[ sqrt(x^2 +y^2) ]

                -> sur (0,y), (L,y), (x,0) et (x,L) avec 0 <= x,y <= L.
                -> sur (x, y)   avec    borninfx <= x <= bornsupx
                                et      borninfy <= y <= bornsupy

    La numérotation des inconnues est lexicographique, la direction x étént
    parcourue avant celle de y. La matrice est retournée dans le format CRS
    qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

    Arguments
    =========
    m           (input)  - nombre de points par direction dans la grille
    n           (output) - pointeur vers le nombre d'inconus dans le système
    L           (input)  - dimensions LxL de la grille (en mètres)
    born_coord  (input) - tableau contenant les coordonées des
                            domaines des parties recatngulaires
    ia          (output) - pointeur vers le tableau 'ia' de la matrice A
    ja          (output) - pointeur vers le tableau 'ja' de la matrice A
    a           (output) - pointeur vers le tableau 'a' de la matrice A
    b           (output) - pointeur vers le tableau 'b'
    nR          (input) - nombre de parties rectangulaires
*/
{
    /* ----- DECLARATION DES VARIABLES -----*/
    int  nnz = 0, ix, iy, nx, ind = 0;
    double h, invh2;
    // nbre de points sur une ligne/colonne de la grille
    // -> noeuds de Dirichlet ne sont pas pris en compte
    nx = m - 2;
    // pas de discrétisation h
    h = L / (m-1);
    // inverse du carré du pas de discrétisation h
    invh2 = 1 / (h*h);
    // nbre d'inconnues
    *n = unknownCalculatorBiss(nx, born_coord, nR, &nnz);

    /* ----- ALOCATION DES TABLEAUX ----- */
    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    // on vérifie si l'allocation est réussite
    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
        return 1;
    }

    /* ----- REMPLISSAGE DE LA MATRICE ----- */
    // initialisation de la première valeur de ia
    nnz = 0;
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            // liste contenant toutes les infos relatives aux
            // conditions aux bords des parties rectangulaires
            int listCond[13] = {0};
            conditionBord(ix, iy, born_coord, nR, nx, listCond);
            
            // on vérifie si on est dans une partie rectangulaire
            if(!listCond[0]) {
                
                // marquer le début de la ligne suivante dans le tableau 'ia'
                (*ia)[ind] = nnz;
                
                // calculer le membre de droite
                (*b)[ind] = 0.0; // rho = 0.0
       
                // remplissage de la ligne : voisin sud
                if ( iy > 0 && !listCond[3] )  {
                    (*a)[nnz] = -invh2; // pour D = 1
                    if(listCond[7]) {
                        (*ja)[nnz] = ind - nx + listCond[8];
                    }
                    else { (*ja)[nnz] = ind - nx; }
                    nnz++;
                }
                else {
                    (*b)[ind] += dirichletCond(ix, iy-1, h) * invh2; // conditions de Dirichlet, bord sud
                }

                // remplissage de la ligne : voisin ouest
                if ( ix > 0 && !listCond[2] )  {
                    (*a)[nnz] = -invh2; // pour D = 1
                    (*ja)[nnz] = ind - 1;
                    nnz++;
                }
                else {
                    (*b)[ind] += dirichletCond(ix-1, iy, h) * invh2; // conditions de Dirichlet, bord ouest
                }
                
                // remplissage de la ligne : élém. diagonal
                (*a)[nnz] = 4.0*invh2; // pour D = 1
                (*ja)[nnz] = ind;
                nnz++;

                // remplissage de la ligne : voisin est
                if ( ix < nx - 1 && !listCond[1] ) {
                    (*a)[nnz] = -invh2; // pour D = 1
                    (*ja)[nnz] = ind + 1;
                    nnz++;
                }
                else {
                    (*b)[ind] += dirichletCond(ix+1, iy, h) * invh2; // conditions de Dirichlet, bord est
                }

                // remplissage de la ligne : voisin nord
                if ( iy < nx - 1 && !listCond[4] ) {
                    (*a)[nnz] = -invh2; // pour D = 1
                    if(listCond[5]) {
                        (*ja)[nnz] = ind + nx - listCond[6];
                    }
                    else { (*ja)[nnz] = ind + nx; }
                    nnz++;
                }
                else {
                    (*b)[ind] += dirichletCond(ix, iy+1, h) * invh2; // conditions de Dirichlet, bord nord
                }
                
                // numéro de l'équation
                ind ++;
                
            }
        }
    }
    
    // dernier élément du tableau 'ia'
    (*ia)[ind] = nnz;    
        
    // retour de fonction habituel
    return 0;
}
