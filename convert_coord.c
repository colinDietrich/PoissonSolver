#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "conditionBord.h"
#define new_max(x,y) (((x) >= (y)) ? (x) : (y))
#define new_min(x,y) (((x) <= (y)) ? (x) : (y))

int convertCoord(double d, double L, int m)
/*
    BUT
    ===
    cette fonction permet de convertir une distance (en  mètre)
    en coordonnée
 
    ARGUMENTS
    =========
    d       (input) - distance à convertir
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                      du domaine de la membrane
    nR      (input) - nombre de parties rectangulaire
    -> renvoie la coordonnée coorespondant à la distance d
 */
{
    int point;
    double h = L/(m-1);
    point = round(d/h - 1);
    return point;
}

void createCoord(int nR, double L, int m, double *born_m, int *born_coord)
/*
    BUT
    ===
    cette fonction permet de convertir un tableau de distance
    (en  mètre) en tableau de coordonnées
    -> plus précisément, cela permet de convertir les positions
    (en mètre) des parties rectangulaires en coordonnées
 
    ARGUMENTS
    =========
    nR      (input) - nombre de parties rectangulaires
    L       (input) - dimensions de la membrane
    m       (input) - rafinement déterminant la discrétisation
                      du domaine de la membrane
    nR      (input) - nombre de parties rectangulaire
    born_m  (input) - tableau contenant les positions (en mètre)
                      des domaines des parties recatngulaires
    born_coord  (output) - tableau contenant les coordonées des
                           domaines des parties recatngulaires
 */
{
    int i = 0;
    for(i = 0; i < nR*4; i++) {
        born_coord[i] = convertCoord(born_m[i], L, m);
    }
}

int calculate_ind(int x, int y, int *born_coord, int nR, int nx)
/*
    BUT
    ===
    cette fonction reçoit les coordonnées (x,y) d'un point sur la grille
    et renvoit l'indice correspondant du vecteur solution du système
    contenant que les points de valeurs non nuls ou non connues.
 
    ARGUMENTS
    =========
    x           (input) - coordonnée selon x du point
    y           (input) - coordonnée selon y du point
    born_coord  (output) - tableau contenant les coordonées des
                           domaines des parties recatngulaires
    nR          (input) - nombre de parties rectangulaire
    nx          (input) - nbre de points sur une ligne/colonne de la grille
 
    => renvoit l'indice de u correspondant au point (x,y)
 */
{
    // incrément
    int i, j, counter = 0;
    int *indexes = malloc((nR + 1) * sizeof(int));
    // nombre de points compris dans les zones rectangulaires
    // entre les coordonnées (0,0) et (ix,iy) sur la fine grid
    int nb = 0;
    
    // on itère sur toutes les zones rectangulaires
    for(i = 0; i < nR; i++) {
        // coordonnées du rectangle
        int x1 = born_coord[4*i], x2 = born_coord[4*i + 1];
        int y1 = born_coord[4*i + 2], y2 = born_coord[4*i + 3];
        
        int bord_x = 0, bord_y = 0;
        
        if(x1 <= 0) { bord_x = 1; }
        if(x2 >= nx) { bord_x = 1; }
        if(y1 <= 0) { bord_y = 1; }
        
        // calcul du nombre de points compris dans les zones rect
        if((x > x2) && (y1 <= y) && (y < y2)) {
            nb += (x2-x1+1 - bord_x) * (y-y1+1 - bord_y);
            indexes[counter+1] = i;
            counter ++;
        }
        else if((x < x1) && (y1 < y) && (y <= y2)) {
            nb += (x2-x1+1 - bord_x) * (y-y1 - bord_y);
            indexes[counter+1] = i;
            counter ++;
        }
        else if( (x>x2 && y==y2) || (y > y2)) {
            nb += (x2-x1+1 - bord_x) * (y2-y1+1 - bord_y);
            indexes[counter+1] = i;
            counter ++;
        }
        indexes[0] = counter;
    }
    
    // on vérfifie que les Ri n'ont pas de cotes en commun
    // si c'est le cas on enleve 1 à nb
    for(i = 0; i < indexes[0]; i++)
    {
        for(j = i+1; j < indexes[0]; j++)
        {
            int xa1 = born_coord[4 * indexes[i+1]], xa2 = born_coord[4 * indexes[i+1] + 1];
            int xb1 = born_coord[4 * indexes[j+1]], xb2 = born_coord[4 * indexes[j+1] + 1];
            int ya1 = born_coord[4 * indexes[i+1] + 2], ya2 = born_coord[4 * indexes[i+1] + 3];
            int yb1 = born_coord[4 * indexes[j+1] + 2], yb2 = born_coord[4 * indexes[j+1] + 3];
            int min = new_min(xa2,xb2);
            int max = new_max(xa1,xb1);
            int minY = new_min(ya2,yb2);
            int maxY = new_max(ya1,yb1);

            // coté haut/bas de Ri avec coté bas/haut de Rj
            if((ya2 == yb1 || ya1 == yb2) && (min-max) > 0)
            {
                int corr1 = (min == nx) ? (-1) : (0);
                int corr2 = (max == -1) ? (-1) : (0);
                nb -= min - max + 1 + corr1 + corr2;
            }
            // coté droit/gauche de Ri avec coté gauche/droit de Rj
            if((xa2 == xb1 || xa1 == xb2) && (minY-maxY) > 0)
            {
                int corr1 = (minY == nx) ? (-1) : (0);
                int corr2 = (maxY == -1) ? (-1) : (0);
                nb -= minY - maxY + 1 + corr1 + corr2;
            }
        }
    }
    int ind = nx * y + x - nb;
    return ind;
}
