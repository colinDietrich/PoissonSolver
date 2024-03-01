#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
double residu(int *ia, int *ja, double *a, double *b, double *u, int n)
/*
    But
    ===
    Génère la norme du résidu du problème : Au = b
    la formule du calcul du résidu est la suivante :
              
                        r_m+1  =  ||b - Au_m+1||
                                  ==============
                                      ||b||

    La matrice A est fournie dans le format CRS qui est défini par le scalaire
    'n' et les trois tableaux 'ia, 'ja' et 'a'.

    Arguments
    =========
    ia          (input) - pointeur vers le tableau 'ia' de la matrice A
    ja          (input) - pointeur vers le tableau 'ja' de la matrice A
    a           (input) - pointeur vers le tableau 'a' de la matrice A
    b           (input) - pointeur vers le tableau 'b' de la
                          matrice colonne b
    u           (input) - pointeur vers l'approximation du tableau
                          du vecteur solution du système Au = b
    n           (input) - nombre d'inconues dans le système
    
    => Renvoit la norme du résidu associée au vecteur u
*/
{
    /* ----- DECLARATION DES VARIABLES ----- */
    // incréments
    int i = 0, j = 0;
    double w = 0.0;
    // square_error = norme carré de l'erreur absolue
    double square_error = 0, square_b = 0, res = 0;
    
    /* ----- CALCUL DU RESIDU ----- */
    for(i = 0; i < n; i++)
    {
        w = 0.0;
        // calcul de Au
        for (j = ia[i]; j < ia[i+1]; j++)
        {
            w += a[j] * u[ ja[j] ];
        }
        // calcul de ||b - Au||^2
        square_error += (b[i] - w) * (b[i] - w);
        // calcul de ||b||^2
        square_b += b[i] * b[i];
    }
    // calcul final du résidu
    res = sqrt(square_error)/sqrt(square_b);
    return res;
}

int residuMatrix(int *ia, int *ja, double *a, double *b, double *u, int n, double *res)
/*
    But
    ===
    Génère la matrice de résidus : r = b - Au

    La matrice A est fournie dans le format CRS qui est défini par le scalaire
    'n' et les trois tableaux 'ia, 'ja' et 'a'.

    Arguments
    =========
    ia          (input) - pointeur vers le tableau 'ia' de la matrice A
    ja          (input) - pointeur vers le tableau 'ja' de la matrice A
    a           (input) - pointeur vers le tableau 'a' de la matrice A
    b           (input) - pointeur vers le tableau 'b' de la
                          matrice colonne b
    u           (input) - pointeur vers l'approximation du tableau
                          du vecteur solution du système Au = b
    n           (input) - nombre d'inconues dans le système
    res         (output) - pointeur vers le tableau 'res' du vecteur résidu
    
    => Renvoit le vecteur résidu associé au vecteur u
*/
{
    /* ----- DECLARATION DES VARIABLES ----- */
    // incréments
    int i = 0, j = 0;
    
    /* ----- CALCUL DU RESIDU ----- */
    for(i = 0; i < n; i++)
    {
        res[i] = 0.0;
        // calcul de Au
        for (j = ia[i]; j < ia[i+1]; j++)
        {
            res[i] -= a[j] * u[ ja[j] ];
        }
        // calcul de b - Au
        res[i] += b[i];
    }
    return 0;
}

double norm(double *u, int n)
/*
    But
    ===
    Calcule la norme du vecteur u

    Arguments
    =========
    u           (input) - pointeur vers le tableau 'u' du vecteur u
    n           (input) - nombre d'éléments du vecteur u
    
    => Renvoit |u|, la norme du vecteur u
*/
{
    /* ----- DECLARATION DES VARIABLES ----- */
    // incréments
    int i = 0;
    double square_u = 0;
    
    /* ----- CALCUL DU RESIDU ----- */
    for(i = 0; i < n; i++)
    {
        // calcul de ||bi||^2
        square_u += u[i] * u[i];
    }
    // calcul final du résidu
    return sqrt(square_u);
}
