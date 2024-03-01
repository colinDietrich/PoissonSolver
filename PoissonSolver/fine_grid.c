#include <stdlib.h>
#include <stdio.h>
#include "conditionBord.h"
#include "dirichlet.h"
#include "write.h"
#include "convert_coord.h"

void prolongation(int m, int n, double L, int *born_coord, double *uR, double *u, int nR)
/*
    But
    ===
    Soit un vecteur 'u' défini sur une grille grossière (coarse-grid) avec un pas
    de discrétisation hR.
    Cette fonction calcule la prolongation 'u' du vecteur 'uR' telle que 'u' soit
    l'approximation du vecteur 'uR' sur une grille plus fine (fine-grid)
    correspondant à un plus petit pas de discrétisation h = hR/2.
 
    La restriction selon les formules suivantes :
                        
        (1) u(2i+1,2j+1)    =       uR(2i+1,2j+1)
        (2) u(2i+1,2j)      = 1/2 [ uR(2i+1,2j-1) + uR(2i+1,2j+1) ]
        (3) u(2i,2j+1)      = 1/2 [ uR(2i-1,2j+1) + uR(2i+1,2j+1) ]
        (4) u(2i,2j)        = 1/4 [ uR(2i-1,2j-1) + uR(2i-1,2j+1)
                                   + uR(2i+1,2j+1) + uR(2i+1,2j-1)]

    Arguments
    =========
    m           (input) - rafinement de la fine-grid
    n           (input) - pointeur vers le nombre d'inconus dans le système
    L           (input) - dimensions du domaine
    born_coord  (input) - pointeur vers le tableau 'born_coord' contenant
                          toutes les coordonnées des parties rectangulaires
                          sur la fine-grid
    uR          (output) - pointeur vers le tableau 'uR' du vecteur réduit uR
    u           (input) - pointeur vers le tableau 'u' du vecteur solution à réduire
    nR          (input) - nombre de parties rectangulaires sur le domaine

    => revoit 0 si fonctionnement normal ou 1 sinon
*/
{
    /* ATTENTION : vecteur de prolongation u doit etre initialiser avec que
     des valeurs nulles */
    
    /* ----- DECLARATION DES VARIABLES -----*/
    int ind = 0, indR = 0, ixR, iyR;
    double w;
    // nbre de points sur une ligne/colonne de la grille
    // -> noeuds de Dirichlet ne sont pas pris en compte
    int nx = m - 2;
    // pas de discrétisation h
    double h = L / (m-1);
    // nbre de points sur une ligne/colonne de la coarse grille
    int nxR;
    if(nx%2 == 0) { nxR = nx/2; }
    else { nxR = (nx-1)/2; }
    
    /* ----- REMPLISSAGE DE VECTEUR U ----- */
    // on itère sur toute la coarse grid
    for (iyR = 0; iyR < nxR; iyR ++) {
        for (ixR = 0; ixR < nxR; ixR++) {
            
            // valeur de ix et iy sur la fine-grid
            int ix = 2*ixR + 1;
            int iy = 2*iyR + 1;
            // liste contenant toutes les infos relatives aux
            // conditions aux bords des parties rectangulaires
            int listCond[13] = {0};
            conditionBord(ix, iy, born_coord, nR, nx, listCond);
            // on vérifie si on est dans une partie rectangulaire
            if(!listCond[0]) {
                w = uR[indR];
                indR++;
                // SUD
                if(!listCond[3]) {
                    ind = calculate_ind(ix, iy-1, born_coord, nR, nx);
                    u[ind] += 0.5*w;
                }
                // SUD-OUEST
                if(!listCond[11]) {
                    ind = calculate_ind(ix-1, iy-1, born_coord, nR, nx);
                    u[ind] += 0.25*w;
                }
                // SUD-EST
                if(!listCond[12]) {
                    ind = calculate_ind(ix+1, iy-1, born_coord, nR, nx);
                    u[ind] += 0.25*w;
                }
                // OUEST
                if(!listCond[2]) {
                    ind = calculate_ind(ix-1, iy, born_coord, nR, nx);
                    u[ind] += 0.5*w;
                }
                // DIAGONALE
                if(!listCond[0]) {
                    ind = calculate_ind(ix, iy, born_coord, nR, nx);
                    u[ind] += w;
                }
                // EST
                if(!listCond[1]) {
                    ind = calculate_ind(ix+1, iy, born_coord, nR, nx);
                    u[ind] += 0.5*w;
                }
                // NORD
                if(!listCond[4]) {
                    ind = calculate_ind(ix, iy+1, born_coord, nR, nx);
                    u[ind] += 0.5*w;
                }
                // NOED-OUEST
                if(!listCond[9]) {
                    ind = calculate_ind(ix-1, iy+1, born_coord, nR, nx);
                    u[ind] += 0.25*w;
                }
                // NORD-EST
                if(!listCond[10]) {
                    ind = calculate_ind(ix+1, iy+1, born_coord, nR, nx);
                    u[ind] += 0.25*w;
                }
            }
        }
    }
}
