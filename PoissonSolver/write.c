#include <stdio.h>

int writeTabInt(int n, int *tab, char name[])
/*
    But
    ===
    Cette fonction permet d'enregistrer une matrice colonne A
    au format int dans un fichier .txt

    Arguments
    =========
    n         (input) - nombre d'éléments du tableau 'tab'
    tab       (input) - tableau de la matrice colonne à enregistrer
    name      (input) - nom du fichier .txt
*/
{
    FILE *file = NULL; int i;
    
    /* ----- ECRITURE "name.txt" ----- */
    file = fopen(name, "w");
    if(file == NULL) {
        printf("\n-----ERREUR-----\n IMPOSSIBLE D'OUVRIR LE FICHIER '%s.txt'", name);
        return 1;
    }
    fprintf(file, "# type: matrix : %s\n# rows : %d\n",name, n);
    for(i=0; i<n; i++) {
        fprintf(file, "%d\n", tab[i]);
    }
    fclose(file);
    
    return 0;
}

int writeTabDouble(int n, double *tab, char name[])
/*
    But
    ===
    Cette fonction permet d'enregistrer une matrice colonne A
    au format int dans un fichier .txt

    Arguments
    =========
    n         (input) - nombre d'éléments du tableau 'tab'
    tab       (input) - tableau de la matrice colonne à enregistrer
    name      (input) - nom du fichier .txt
*/
{
    FILE *file = NULL; int i;
    
    /* ----- ECRITURE "name.txt" ----- */
    file = fopen(name, "w");
    if(file == NULL) {
        printf("\n-----ERREUR-----\n IMPOSSIBLE D'OUVRIR LE FICHIER '%s.txt'", name);
        return 1;
    }
    fprintf(file, "# type: matrix : %s\n# rows : %d\n",name, n);
    for(i=0; i<n; i++) {
        fprintf(file, "%e\n", tab[i]);
    }
    fclose(file);
    
    return 0;
}

int appendText(double value, char name[])
/*
    But
    ===
    Cette fonction permet d'enregistrer une matrice colonne A
    au format int dans un fichier .txt

    Arguments
    =========
    n         (input) - nombre d'éléments du tableau 'tab'
    tab       (input) - tableau de la matrice colonne à enregistrer
    name      (input) - nom du fichier .txt
*/
{
    FILE *file = NULL; int i;
    
    /* ----- ECRITURE "name.txt" ----- */
    file = fopen(name, "a");
    if(file == NULL) {
        printf("\n-----ERREUR-----\n IMPOSSIBLE D'OUVRIR LE FICHIER '%s.txt'", name);
        return 1;
    }
    fprintf(file, "%e\n", value);
    fclose(file);
    return 0;
}

int writeCSR_to_Matrix(int n, int nia, int *ia, int *ja, double *a, char name[])
{
    FILE *file = NULL;
    
    file = fopen(name, "w");
    if(file == NULL) {
        printf("\n-----ERREUR-----\n IMPOSSIBLE D'OUVRIR LE FICHIER '%s.txt'", name);
        return 1;
    }
    fprintf(file, "# type: matrix : %s\n# rows : %d\n# columns : %d\n", name, n, n);
    
    // incréments d'itération
    int i, j, col, inc = 0;
    
    for(i = 0; i < n; i++) {
        // numéro de colonne
        col = 0;
        for(j = 0; j < n; j++) {
            if(j == ja[ia[i]+col]) {
                fprintf(file, "%10.2f", a[ia[i]+col]);
                if(col < ia[i+1]-ia[i]-1) {
                    col++;
                }
            } else {
                fprintf(file, "%10.2f",0.0);
            }
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    
    return 0;
}
