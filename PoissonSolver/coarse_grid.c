#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "write.h"
#include "unknown_calculator.h"
#include "conditionBord.h"
#include "dirichlet.h"
#include "convert_coord.h"

void restriction (int m, int n_R, int *born_coord, int *born_coordR, int nR, double *uR, double *u, double L)
/*
    But
    ===
    Soit un résidu 'u' défini sur une grille fine (fine-grid) avec un pas
    de discrétisation h.
    Cette fonction calcule la restriction 'uR' du vecteur 'u' telle que 'uR' soit
    l'approximation du vecteur 'u' sur une grille plus grossière (coarse-grid)
    correspondant à un plus grand pas de discrétisation h_r = 2*h.
 
    La restriction esr réalisée par FULL-WEIGHTING :
                        
    V(2i+1,2j+1) =    1/4  *   u(2i+1,2j+1)
                    + 1/8  * [ u(2i,2j+1) + u(2i+2,2j+1) + u(2i+1,2j) + u(2i+1,2j+2) ]
                    + 1/16 * [ u(2i,2j)   + u(2i+2,2j)   + u(2i,2j+2) + u(2i+2,2j+2) ]

    Arguments
    =========
    m           (input) - rafinement de la fine-grid
    n_R         (input) - nombre d'éléments du vecteur réduit 'uR'
    born_coord  (input) - pointeur vers le tableau 'born_coord' contenant
                          toutes les coordonnées des parties rectangulaires
                          sur la fine-grid
    born_coordR (input) - pointeur vers le tableau 'born_coordT' contenant
                          toutes les coordonnées des parties rectangulaires
                          sur la coarse-grid
    nR          (input) - nombre de parties rectangulaires sur le domaine
    uR          (output) - pointeur vers le tableau 'uR' du vecteur réduit uR
    u           (input) - pointeur vers le tableau 'u' du vecteur solution à réduire
    L           (input) - dimensions du domaine

    => revoit 0 si fonctionnement normal ou 1 sinon
*/
{

    /* ----- DECLARATION DES VARIABLES -----*/
    // indices du vecteur u
    int ind = 0;    // fine grid
    int indR = 0;   // coarse grid
    // indices (x,y)
    int ix, iy;     // fine grid
    int ixR, iyR;   // coarse grid
    // nbre de points sur une ligne/colonne de la grille
    int nx = m - 2; // fine grid
    int nxR;        // coarse grid
    if(nx%2 == 0) { nxR = nx/2; }
    else { nxR = (nx-1)/2; }
    // pas de discrétisation h
    double h = L / (m-1);
    // mise en mémoire de la composante i du vecteur uR
    double w;
    // poids des différents points voisins
    double c1 = 1.0/4.0, c2 = 1.0/8.0, c3 = 1.0/16.0;
    int index;
    
    /* ----- ITERATION SUR LA COARSE-GRID ----- */
    for (iyR = 0; iyR < nxR; iyR ++) {
        for (ixR = 0; ixR < nxR; ixR++) {
            ix = 2*ixR + 1;
            iy = 2*iyR + 1;
            // liste contenant toutes les infos relatives aux
            // conditions aux bords des parties rectangulaires
            int listCond[13] = {0};
            conditionBord(ix, iy, born_coord, nR, nx, listCond);
            
            // on vérifie si on est dans une partie rectangulaire
            if(!listCond[0]) {
                ind = calculate_ind(ix, iy, born_coord, nR, nx);
                // réinitialisation de la variable de stockage
                w = 0.0;
                
                // calcul de l'indice du point sud sur la fine-grid
                if(listCond[7]) {
                    index = ind - nx + listCond[8];
                } else {
                    index = ind - nx;
                }
                // SUD
                if(iy > 0 && !listCond[3]) {
                    w += c2 * u[index];
                }
                // SUD-OUEST
                if (ix-1 >= 0 && iy > 0 && !listCond[11]) {
                    w += c3 * u[index - 1];
                }
                // SUD-EST
                if (ix+1 < nx && iy > 0 && !listCond[12]) {
                    w += c3 * u[index + 1];
                }
                // OUEST
                if ( ix > 0 && !listCond[2] )  {
                    w += c2 * u[ind - 1];
                }
                // DIAGONALE
                w += c1 * u[ind];
                // EST
                if ( ix < nx - 1 && !listCond[1] ) {
                    w += c2 * u[ind + 1];
                }
                // calcul de l'indice du point nord sur la fine-grid
                if(listCond[5]) {
                    index = ind + nx - listCond[6];
                }
                // NORD
                if(iy < nx - 1 && !listCond[4]) {
                    w += c2 * u[index];
                }
                // NORD-OUEST
                if (ix-1 >= 0 && iy < nx - 1 && !listCond[9]) {
                    w += c3 * u[index - 1];
                }
                // NORD-EST
                if (ix+1 < nx && iy < nx-1 && !listCond[10]) {
                    w += c3 * u[index + 1];
                }

                // mise en mémoire de la composante i du vecteur réduit
                uR[indR] = w;
                
                // numéro de l'équation
                indR ++;
            }
        }
    }
}
