#include <stdio.h>
#include <math.h>

void add(double *u, double *corr, int n, int niv)
/*
    BUT
    ===
    cette fonction somme le vecteur u au vecteur corr si l'entier
    niv est différent de 0 et fait l'inverse le cas échéant
 
    ARGUMENTS
    =========
    u       (input) - vecteur à sommer
    corr    (input) - vecteur à sommer
    n       (input) - longueur des vecteurs u/corr
    niv     (input) - condition de sommation
*/
{
    int i;
    if(niv != 0) {
        for(i = 0; i < n; i++) {
            corr[i] += u[i];
        }
        
    } else {
        for(i = 0; i < n; i++) {
            u[i] += corr[i];
        }
    }
}

void clearArray(double *v, int n)
/*
    BUT
    ===
    cette fonction réinitialise un vecteur en mettant
    toutes ses valeurs sur 0
 
    ARGUMENTS
    =========
    v       (input) - vecteur à sréinitialiser
    n       (input) - longueur du vecteur v
*/
{
    int i;
    for(i = 0; i < n; i++) {
        v[i] = 0.0;
    }
}

void update(double *u, double *uOld, int n, double tau)
/*
    BUT
    ===
    cette fonction met à jour le vecteur solutiion calculer par
    l'algorithme Multi-grid en prenant en compte le facteur de
    relaxation
 
    ARGUMENTS
    =========
    u          (input) - vecteur solution à mettre à jour
    uOld       (input) - vecteur contenant les anciennes valeurs de u
    n          (input) - longeur des vecteur u/uOld
    tau        (input) - facteur de relaxation
*/
{
    int i;
    for(i = 0; i < n; i++) {
        u[i] = uOld[i] + tau * (u[i] - uOld[i]);
    }
}

void copy(double *u, double *u2, int n)
/*
    BUT
    ===
    cette fonction copy un vecteur dans un autre
 
    ARGUMENTS
    =========
    u          (input) - vecteur contenant la copie
    u2         (input) - vecteur à copier
    n          (input) - longeur des vecteur u/u2
*/
{
    int i;
    for(i = 0; i < n; i++) {
        u[i] = u2[i];
    }
}
