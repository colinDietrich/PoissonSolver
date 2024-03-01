#include <stdio.h>
#include <stdlib.h>
#define new_max(x,y) (((x) >= (y)) ? (x) : (y))
#define new_min(x,y) (((x) <= (y)) ? (x) : (y))

int condR(int x, int y, int *born_coord, int nR){
    /*
        BUT
        ===
        cette fonction permet de déterminer si un point de
        coordonnées (x,y) est compris dans une partie rectangulaire
     
        ARGUMENTS
        =========
        x       (input) - abscisse du point
        y       (input) - oordonnée du point
        born_coord  (input) - tableau contenant les coordonées des
                              domaines des parties recatngulaires
        nR      (input) - nombre de parties rectangulaire
        -> renvoie 1 si le point est contenu dans une partie
           rectangulaire, 0 sinon.
     */
    int i = 0, bool = 0;
    for(i = 0; i < nR; i++) {
        if ( (x>=born_coord[4*i] && x<=born_coord[4*i + 1]) && (y>=born_coord[4*i + 2] && y<=born_coord[4*i + 3]) ) {
            bool = 1;
        }
    }
    return bool;
}

int nbPointsInR(int *born_coord, int nx, int *indexes)
/*
    BUT
    ===
    cette fonction permet de déterminer le nombre de points
    contenu dans le Ri par lesquels il faut passer pour atteindre
    le voisin nord/sud
 
    ARGUMENTS
    =========
    born_coord  (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires
    nx          (input) - entier = m-2
    indexes     (input) - tableau contenant le nombre et les
                          indices des parties rect Ri précisées
                          ci-dessus
 */
{
    int i = 0, j = 0, nb = 0;
    for(i = 0; i < indexes[0]; i++)
    {
        // nbre de point contenu dans Ri
        nb += born_coord[4 * indexes[i+1] + 1] - born_coord[4 * indexes[i+1]] + 1;
        
        // si un x1 de Ri touche le bord de la membrane carrée initiale
        // -> on enleve 1 a nb
        if(born_coord[4 * indexes[i+1]] == -1 || born_coord[4 * indexes[i+1]] == nx)
        {
            nb -= 1;
        }
        // si un x2 de Ri touche le bord de la membrane carrée initiale
        // -> on enleve 1 a nb
        if(born_coord[4 * indexes[i+1] + 1] == -1 || born_coord[4 * indexes[i+1] + 1] == nx)
        {
            nb -= 1;
        }
        
        // on vérfifie que les Ri n'ont pas de cotes en commun
        // si c'est le cas on enleve 1 à nb
        for(j = i+1; j < indexes[0]; j++)
        {
            int xa1 = born_coord[4 * indexes[i+1]], xa2 = born_coord[4 * indexes[i+1] + 1];
            int xb1 = born_coord[4 * indexes[j+1]], xb2 = born_coord[4 * indexes[j+1] + 1];
            int min = new_min(xa2,xb2);
            int max = new_max(xa1,xb1);
            // coté droit/gauche de Ri avec coté gauche/droit de Rj
            if(xa2 == xb1 || xa1 == xb2)
            {
                nb -= 1;
            }
            // coté haut/bas de Ri avec coté bas/haut de Rj
            if((born_coord[4 * indexes[i+1] + 3] == born_coord[4 * indexes[j+1] + 2] || born_coord[4 * indexes[i+1] + 2] == born_coord[4 * indexes[j+1] + 3]) && (min-max) > 0)
            {
                int corr1 = (min == nx) ? (-1) : (0);
                int corr2 = (max == -1) ? (-1) : (0);
                nb -= min - max + 1 + corr1 + corr2;
            }
        }
    }
    return nb;
}


void conditionBord(int x, int y, int *born_coord, int nR, int nx, int listCond[13])
/*
    BUT
    ===
    cette fonction permet de déterminer les conditions suivantes
    pour un point de coordonnées (x,y) :
    (1) listCond[0] = 1 si le point est compris dans une partie
        rectangulaire, = 0 sinon
    (2) listCond[1] = 1 si le point est situé à côté du bord
        gauche (ouest) d'une partie rectangulaire, = 0 sinon
    (3) listCond[2] = 1 si le point est situé à côté du bord
        de droite (est) d'une partie rectangulaire, = 0 sinon
    (4) listCond[3] = 1 si le point est situé au dessus du bord
        du haut (nord) d'une partie rectangulaire, = 0 sinon
    (5) listCond[4] = 1 si le point est situé en bas du bord du
        bas (sud) d'une partie rectangulaire, = 0 sinon
    (6) listCond[5] = 1 si le chemin à parcourirs pour aller
        au voisin nord passe par des parties rectangulaires,
        = 0 sinon
    (7) listCond[6] = le nombre de points contenu dans le Ri par
        lesquels il faut passer pour atteindre le voisin nord
    (8) listCond[7] = 1 si le chemin à parcourirs pour aller au
        voisin sud passe par des parties rectangulaires
    (9) listCond[8] = le nombre de points contenu dans le Ri par
        lesquels il faut passer pour atteindre le voisin sud
    (10) listCond[9] = 1 si voisin nord-ouest existe
    (11) listCond[10] = 1 si voisin nord-est existe
    (12) listCond[11] = 1 si voisin sud-ouest existe
    (13) listCond[12] = 1 si voisin sud-est existe
 
    ARGUMENTS
    =========
    x       (input) - abscisse du point
    y       (input) - oordonnée du point
    born_coord  (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires
    nR      (input) - nombre de parties rectangulaire
    nx      (input) - entier = m-2
    listCond    (output) - voir ci-dessus
 */
{
    int i = 0, counterSouth = 0, counterNorth = 0;
    // détermination du nombre de points contenu dans les Ri
    // par lesquels il faut passer pour atteindre le voisin sud/nord
    int *indexesSouth = malloc((nR + 1) * sizeof(int));
    int *indexesNorth = malloc((nR + 1) * sizeof(int));
    
    // on vérifie les conditions aux bords sur tous les Ri
    for(i = 0; i < nR; i++) {
        
        int x1 = born_coord[4*i], x2 = born_coord[4*i + 1];
        int y1 = born_coord[4*i + 2], y2 = born_coord[4*i + 3];
        
        // déterminer si un point de coordonnées (x,y) est compris dans
        // une partie rectangulaire
        if ( (x>=x1 && x<=x2) && (y>=y1 && y<=y2) ) {
            listCond[0] = 1;
        }
        // déterminer si un point de coordonnées (x,y) est situé à
        // côté du bord gauche (ouest) d'une partie rectangulaire
        if( (x+1>=x1 && x+1<=x2) && (y>=y1 && y<=y2) ) {
            listCond[1] = 1;
        }
        // déterminer si un point de coordonnées (x,y) est situé à
        // côté du bord de droite (est) d'une partie rectangulaire
        if( (x-1>=x1 && x-1<=x2) && (y>=y1 && y<=y2) ) {
            listCond[2] = 1;
        }
        // déterminer si un point de coordonnées (x,y) est situé au
        // dessus du bord du haut (nord) d'une partie rectangulaire
        if( (x>=x1 && x<=x2) && (y-1>=y1 && y-1<=y2) ) {
            listCond[3] = 1;
        }
        // déterminer si un point de coordonnées (x,y) est situé en
        // bas du bord du bas (sud) d'une partie rectangulaire
        if( (x>=x1 && x<=x2) && (y+1>=y1 && y+1<=y2) ) {
            listCond[4] = 1;
        }
        // détermine si le chemin à parcourirs pour aller au voisin nord
        // passe par des parties rectangulaires + enregistre le nombres
        // et les indices des parties rectangulaires
        if( ((y>=y1-1 && y<y2) && (x>=x2 && x<=nx-1)) || ((y>=y1 && y<=y2) && (x>=0 && x<=x1)) ) {
            listCond[5] = 1;
            indexesNorth[counterNorth + 1] = i;
            counterNorth ++;
        }
        indexesNorth[0] = counterNorth;
        // détermine si le chemin à parcourirs pour aller au voisin sud
        // passe par des parties rectangulaires + enregistre le nombres
        // et les indices des parties rectangulaires
        if( ((y>=y1 && y<=y2) && (x>=x2 && x<=nx-1)) || ((y>=y1+1 && y<=y2+1) && (x>=0 && x<=x1)) )
        {
            listCond[7] = 1;
            indexesSouth[counterSouth + 1] = i;
            counterSouth ++;
        }
        indexesSouth[0] = counterSouth;
        
        // détermine si voisin nord-ouest existe
        if ((x-1>=x1 && x-1<=x2) && (y+1>=y1 && y+1<=y2)) {
            listCond[9] = 1;
        }
        // détermine si voisin nord-est existe
        if ((x+1>=x1 && x+1<=x2) && (y+1>=y1 && y+1<=y2)) {
            listCond[10] = 1;
        }
        // détermine si voisin sud-ouest existe
        if ((x-1>=x1 && x-1<=x2) && (y-1>=y1 && y-1<=y2)) {
            listCond[11] = 1;
        }
        // détermine si voisin sud-est existe
        if ((x+1>=x1 && x+1<=x2) && (y-1>=y1 && y-1<=y2)) {
            listCond[12] = 1;
        }
    }
    // détermine le nombre de points contenu dans le Ri par
    // lesquels il faut passer pour atteindre le voisin nord
    listCond[6] = nbPointsInR(born_coord, nx, indexesNorth);
    // détermine le nombre de points contenu dans le Ri par
    // lesquels il faut passer pour atteindre le voisin sud
    listCond[8] = nbPointsInR(born_coord, nx, indexesSouth);
    
    // libérer la mémoire
    free(indexesNorth); free(indexesSouth);
}
