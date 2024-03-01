
#include <stdio.h>
#include "conditionBord.h"



int unknownCalculator(int nx, int *born_coord, int nR)
/*
    But
    ===
    (1) cette fonction permet de déterminer le nombre d'inconnues
        du sytème : Ax = b
 
    (2) n = omega_i - rect + omega_rect
        avec
        -> omega_i = (nx*nx) = nombre de points contenu dans la
           membrane carrée initiale de dimension (LxL)
        -> rect = nombre de points contenus dans les domaines
           de toutes les parties rectangulaires
        -> omega_rect = nombre de points à l'intersection
           d'une partie rectangulaire et du bord de la
           membrane carrée initiale de dimension (LxL)
 
    !!! ATTENTION !!!
    cette fonction ne marche que pour nR = 1

    Arguments
    =========
    nx          (input)  - nombre de points par direction dans la
                           grille sans les bord : nx = m-2
    born_coord  (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires
    nR          (input) - nombre de parties rectangulaires
 
    -> Renvoit le nombre d'inconnues du système
*/
{
    // incrément
    int i = 0;
    // nombre de points contenus dans les domaines
    // de toutes les parties rectangulaires
    int rect = 0;
    // nombre de points contenu dans la
    // membrane carrée initiale de dimension (LxL)
    int omega_i = nx * nx;
    // nombre de points à l'intersection
    // d'une partie rectangulaire et du bord de la
    // membrane carrée initiale de dimension (LxL)
    int omega_rect = 0;
    // nombre d'inconnues du système
    int n = 0;
    // booléen
    int condx1, condx2, condy1, condy2;
    // coordonnées des bords du domaine de la partie rectangulaire Ri
    int x1, x2, y1, y2;
    
    for(i = 0; i < nR; i++) {
        
        x1 = born_coord[4*i];
        x2 = born_coord[4*i + 1];
        y1 = born_coord[4*i + 2];
        y2 = born_coord[4*i + 3];
        
        // si un x1 de Ri touche le bord de la membrane carrée initiale
        condx1 = 0;
        if(x1 == -1 || x1 == nx) { condx1 = 1; }
        // si un x2 de Ri touche le bord de la membrane carrée initiale
        condx2 = 0;
        if(x2 == -1 || x2 == nx) { condx2 = 1; }
        // si un y1 de Ri touche le bord de la membrane carrée initiale
        condy1 = 0;
        if(y1 == -1 || y1 == nx) { condy1 = 1; }
        // si un y2 de Ri touche le bord de la membrane carrée initiale
        condy2 = 0;
        if(y2 == -1 || y2 == nx) { condy2 = 1; }

        
        // on calcule le nombre de points contenu dans chacun des Ri
        rect += (x2 - x1 + 1) * (y2 - y1 + 1);
        
        if( (condx1 && condx2 && condy1) || (condx1 && condx2 && condy2))
        {
            omega_rect += 2 * (y2 - y1) + x2 - x1 + 1;
        }
        else if( (condy1 && condy2 && condx1) || (condy1 && condy2 && condx2) )
        {
            omega_rect += 2 * (x2 - x1) + y2 - y1 + 1;
        }
        else if( (condx1 || condx2) && (condy1 || condy2) )
        {
            omega_rect += y2 - y1 + x2 - x1 + 1;
        }
        else if(condx1 && condx2)
        {
            omega_rect += 2 * (y2 - y1 + 1);
        }
        else if(condy1 && condy2)
        {
            omega_rect += 2 * (x2 - x1 + 1);
        }
        else if(condx1 || condx2)
        {
            omega_rect += y2 - y1 + 1;
        }
        else if(condy1 || condy2)
        {
            omega_rect += x2 - x1 + 1;
        }
        
    }

    // calcul final du nombre d'inconnues du sytème
    n = omega_i - rect + omega_rect;
    return n;
}


int unknownCalculatorBiss(int nx, int *born_coord, int nR, int *nnz)
/*
    But
    ===
    (1) cette fonction permet de déterminer le nombre d'inconnues
        du sytème : Ax = b
    
    !!! ATENTION !!!
    cette fonction, bien que moins efficace que 'unknownCalculator(..)',
    fonctionne pour nR >= 1

    Arguments
    =========
    nx          (input)  - nombre de points par direction dans la
                           grille sans les bord : nx = m-2
    born_coord  (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires
    nR          (input) - nombre de parties rectangulaires
 
    -> Renvoit le nombre d'inconnues du système
*/
{
    int i = 0, ix = 0, iy = 0, n = 0;
    for(i = 0; i < nx*nx; i++) {
        if(ix == nx) {
            iy++;
            ix = 0;
        }
        
        // liste contenant toutes les infos relatives aux
        // conditions aux bords des parties rectangulaires
        int listCond[13] = {0};
        conditionBord(ix, iy, born_coord, nR, nx, listCond);
        
        // on vérifie si on est dans une partie rectangulaire
        if(!listCond[0]) {
            n++;
            // voisin sud
            if ( iy > 0 && !listCond[3] )  { (*nnz)++; }
            // voisin ouest
            if ( ix > 0 && !listCond[2] )  { (*nnz)++; }
            // élém. diagonal
            (*nnz)++;
            // voisin est
            if ( ix < nx - 1 && !listCond[1] ) { (*nnz)++; }
            // voisin nord
            if ( iy < nx - 1 && !listCond[4] ) { (*nnz)++; }
        }
        ix++;
    }
    return n;
}

