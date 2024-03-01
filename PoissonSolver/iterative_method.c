#include <stdlib.h>
#include <stdio.h>
#include "residu.h"

void gaussSeidel(int *ia, int *ja, double *a, double *b, int n, double *u, int upper)
/*
    But
    ===
    Calcule le vecteur solution du système linéaire (n x n)
                                Au = b
    par Par la méthode itérative de Gauss-Seidel inférieure / supérieur

    Le préconditionneur correspondant est la matrice triangilaire
    inférieure / supérieur
                                B = L       si upper = 0
                                B = U       si upper = 1
 
    Algorithme sous forme matricielle :
             (1) Initialisation    -> u_0 = 0 (vecteur nul)
                                   -> r_0 = b - Au_0

             (2) For m = 0,1,...,iteration :
                                  -> d_m      =   B^-1 * r_m
                                  -> u_m+1    =   u_m + d_m
                                  -> r_m+1    =   b - Au_m+1

                => u_k+1 = t.B^-1 (b - Au_k)

    Arguments
    =========
    ia          (inupt) - pointeur vers le tableau 'ia' de la matrice A
    ja          (inupt) - pointeur vers le tableau 'ja' de la matrice A
    a           (inupt) - pointeur vers le tableau 'a' de la matrice A
    b           (inupt) - pointeur vers le tableau 'b'
    n           (input) - pointeur vers le nombre d'inconus dans le système
    u           (output) - pointeur vers le tableau 'u' du vecteur solution
    upper       (input) - booléen -> si 0 : B = L
                                  -> si 1 : B = U
*/
{
    // incréments
    int iter, i, j;
    // On définit les paramètres pour GS INFERIEUR par défaut
    int sign = 1;
    int init = 0;
    int end = n;
    double w, diag;
    // si upper = 1 --> On définit les paramètres pour GS SUPERIEUR
    if(upper) {
        sign = -1;
        init = n-1;
        end = 1;
    }
    
    i = 0, j = 0;
    for(i = init; sign*i < end; i += 1*sign) {
        w = b[i];
        for (j = ia[i]; j < ia[i+1]; j++)
        {
            if(ja[j] != i) {
                w -= a[j] * u[ ja[j] ];
            }
            else {
                diag = a[j];
            }
        }
        w /= diag;
        u[i] = w; // on met à jour la composante i de u_k+1
    }
}


void gaussSeidelSymetrique(int *ia, int *ja, double *a, double *b, int n, double *u)
/*
    But
    ===
    Cette fonction implémente la méthode itérative de Gauss-Seidel symétrique
    à partir de la fonction 'gaussSeidel' en alternant une méthode itérative
    supérieur puis inférieure

    Arguments
    =========
    ia          (inupt) - pointeur vers le tableau 'ia' de la matrice A
    ja          (inupt) - pointeur vers le tableau 'ja' de la matrice A
    a           (inupt) - pointeur vers le tableau 'a' de la matrice A
    b           (inupt) - pointeur vers le tableau 'b'
    n           (input) - pointeur vers le nombre d'inconus dans le système
    u           (output) - pointeur vers le tableau 'u' du vecteur solution
*/
{
    // Gauss-Seidel supérieur
    gaussSeidel(ia, ja, a, b, n, u, 1);
    // Gauss-Seidel inférieur
    gaussSeidel(ia, ja, a, b, n, u, 0);
    
}
