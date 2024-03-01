#include <math.h>
#include <stdio.h>

double dirichletCond(int ix, int iy, double h)
{
/*
    Cette fonction prend en entrée les coordonnées d'un point
    du bord du domaine et renvoie les conditions aux bords de
    Dirichlet :     u = exp[ sqrt(x^2 +y^2) ]
    
     Arguments
     =========
     ix         (input)  - coordonnée en x du point
     iy         (input) - coordonnée en y du point
     h          (input)  - pas de discrétisation du domaine
    
    => renvoit la condition de dirichlet correspondant au point
        situé sur le bord du domaine
 */
    // calcul x^2 et y^2
    double xSquared = (((double)ix+1.0) * h) * (((double)ix+1.0) * h);
    double ySquared = (((double)iy+1.0) * h) * (((double)iy+1.0) * h);
    // renvoie condition de dirichlet sur le bord du domaine
    return exp(sqrt(xSquared + ySquared));
}
