#include <stdio.h>
#include <math.h>
#include "conditionBord.h"

int printVector(double *evecs, double L, int *born_coord, int m, int n, int nR, char name[])
/*
    BUT
    ===
    cette fonction permet d'enregistrer le vecteur u(x,y)
    solution du système dans un fichier texte.
    + afficher les valeurs de u(x,y) avec GNUPlot
 
    ARGUMENTS
    =========
    evecs       (input) - le tableau des vecteurs propres de la
                          membrane (éléments du vecteurs solution
                          u(x,y) non nuls)
    L           (input) - dimensions LxL de la grille (en mètres)
    born_coord  (input) - tableau contenant les coordonées des
                          domaines des parties recatngulaires
    m           (input)  - nombre de points par direction dans
                           la grille
    n           (output) - pointeur vers le nombre d'inconus dans le
                           système
    nR          (input) - nombres de parties rectangulaires
 
    -> Retourne 0 si l'opération s'est bien déroulé, 1 si non.
 */
{
    
    int ind = 0, iy = 0, ix = 0;;
    // pas de discrétisation h
    double h = L / (m-1);
    
    // fichier contenant les valeurs du vecteur u(x,y)
    // solution du système
    FILE *file = NULL;
    file = fopen(name, "w");
    if(file == NULL) {
        printf("\n-----ERREUR-----\n IMPOSSIBLE D'OUVRIR LE FICHIER 'oscillation.txt'");
        return 1;
    }
    
    /* ----- ECRITURE DU VECTEUR DANS FICHIER TEXTE ----- */
    fprintf(file, "#Y\t\tX\t\tUij\n");
    for (iy = 0; iy < m; iy++)
    {
        for (ix = 0; ix < m; ix++)
        {
            // sin on est sur le bord de la membrane -> Uij = 0
            if(ix == 0 || ix == m-1 || iy == 0 || iy == m-1){
                fprintf(file, "%e\t%e\t%e\n", iy*h, ix*h, 0.0);
            }
            // si on est ni dans la zone rectangulaire ni sur son bord
            // -> Uij != 0
            else if(!condR(ix-1, iy-1, born_coord, nR))
            {
                fprintf(file, "%e\t%e\t%e\n", iy*h, ix*h, fabs(evecs[ind]));
                //fprintf(file, "%e\t%e\t%e\n", iy*h, ix*h, 2.0);
                ind++;
            }
            // si on est dans la zone rectangulaire bord inclu -> Uij = 0
            else {
                fprintf(file, "%e\t%e\t%e\n", iy*h, ix*h, 0.0);
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
    
    /* ----- CREATION DU FICHIER DE COMMANDE POUR GNUPLOT ----- */
    /*
    //Création du fichier de commandes
    FILE *gnuplot1 = popen("gnuplot -persistent", "w");
    if(gnuplot1 == NULL) {
        printf("\n-----ERREUR-----\nIMPOSSIBLE D'OUVRIR LE FICHIER DANS LA CREATION DU FICHIER DE COMMANDE POUR GNUPLOT\n");
        return 1;
    }
    
    fprintf(gnuplot1, "load 'cmdGnuplot.txt'\n");
    //fprintf(gnuplot1, "set pm3d map\n");
    fprintf(gnuplot1, "set xrange[%f:%f]\n", 0.0, L);
    fprintf(gnuplot1, "set yrange[%f:%f]\n", 0.0, L);
    fprintf(gnuplot1, "splot 'vec.txt' using 2:1:3 with pm3d\npause mouse close\n");
    fclose(gnuplot1);
     */
     
    return 0;
}
