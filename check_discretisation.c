#include <stdio.h>
#include <math.h>

int checkDiscretisation(int m, double L, double *born_m, int nR, double tol)
/*
    BUT
    ===
    cette fonction détermine si le rafinement m permet de
    discrétiser correctement le domaine ou non.
    Pour cela, il faut que la condition suivante soit vérifiée :
 
    -> bord_m % h < tol
 
    avec -> bord_m = le tableau contenant les cordonnées des bords
         des parties rectangulaires
         -> h = L / (m-1) le pas de discrétisation du domaine
 
    ARGUMENTS
    =========
    m       (input) - entier correspondant au rafinement de la discrétisation
    L       (input) - dimensions de la membrane
    born_m  (input) - tableau contenant les domaines des parties
                      recatngulaires
    nR      (input) - nombre de parties rectangulaire
    tol     (input) - la tolérance déterminant la précision des points
                      de discrétisation pour les bords de la membrane
 
    -> Renvoit 0 si la discrétisation est bien déterminée, 1 sinon.
 */
{
    // incrément
    int i = 0;
    // on vérifie que chaque bord du domaine des parties rectangulaires
    // peut être défini sous une certaine tolérance
    for(i = 0; i < nR*4; i++) {
        
        double rest = born_m[i] / (L / (double)(m-1));
        
        //printf("\n-- born_m[%d] : %e", i, born_m[i]);
        //printf("\nrest : %e", rest);
        //printf("\nround : %e", round(rest));

        // si l'écart de discrétisation par rapport au bord
        // est plus grand que tol -> on empêche l'éxécution
        if( (round(rest)-rest > tol) || (round(rest)-rest < -tol) ) {
            printf("\n\n----- ERREUR : LA VALEUR DU PARAMETRE M NE PERMET PAS DE DISCRETISER LE DOMAINE CORRECTEMENT -----\n\n");
            return 1;
        }
    }
    return 0;
}

int checkMultigrid(int m, double L, double *born_m, int nR, double tol, int levels)
{
/*
    BUT
    ===
    cette fonction détermine si le nombre de niveux de la méthode
    multgrid permet de discrétiser correctement le domaine ou non.
    Pour ce fiare, on calcule le rafinement m de la coarsest grid
    et on utilise ensuite la fonction 'checkDiscretisation' pour
    vérifier que la discrétisation correspondante est autorisée ou non.
 
    ARGUMENTS
    =========
    m       (input) - entier correspondant au rafinement de la discrétisation
    L       (input) - dimensions de la membrane
    born_m  (input) - tableau contenant les domaines des parties
                      recatngulaires
    nR      (input) - nombre de parties rectangulaire
    tol     (input) - la tolérance déterminant la précision des points
                      de discrétisation pour les bords de la membrane
    levels  (input) - nombre de niveux de l'algorithme multigrid
 
    -> Renvoit 0 si la discrétisation est bien déterminée, 1 sinon.
*/
    
    // discrétisation de la coarsest grid
    int final_m;
    int power = (int)pow(2.0, (double)(levels-1));
    if((m + power - 1) % power != 0) {
        printf("\n\n----- LE NOMBRE DE NIVEAU EST TROP ELEVE POUR DISCRETISER COORECTEMENT LE DOMAINE A LA COARSEST-GRID -----\n\n");
        return 1;
    } else {
        final_m = (m + pow(2, levels-1) - 1) / pow(2, levels-1);
    }
    
    // on vérifie que la coarsest grid permet de définir toutes les zones rect
    if(checkDiscretisation(final_m, L, born_m, nR, tol)) {
        printf("\n\n----- LE NOMBRE DE NIVEAU EST TROP ELEVE POUR DISCRETISER COORECTEMENT LE DOMAINE A LA COARSEST-GRID -----\n\n");
        return 1;
    }
    return 0;
}
