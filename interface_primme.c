#include <stdlib.h>
#include <math.h>
#include "primme.h"
#include "check_discretisation.h"
#include "multiGrid.h"

/* variables statiques -- accessibles  */
static double L, *a, *born_m;
static int n, m, nR, *ia, *ja;

void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme)
/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx
   pour le solveur aux valeurs propres PRIMME. La matrice A doit être
   "stoquée" au préalable dans les variables statiques 'n', 'ia', 'ja' et 'a'
   en utilisant le format CSR (Compressed Sparse Rows). Par "stoquer"
   on veut dire ici stoquer la valeur de 'n' et les pointeurs vers les
   tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   vy       (output) - vecteur(s) de produit A*vx
   blockSize (input) - nombre de vecteurs d'entrée
   primme    (input) - paramètres fournis par primme pour optimiser le calcul
                       (pas utilisé)
*/
{
    int i, j, b;
    double *x = vx, *y=vy;

    for(b = 0; b < (*blockSize)*n; b+=n)
        for(i = 0; i < n; i++){
            y[b+i] = 0;
            for (j = ia[i]; j < ia[i + 1]; j++)
                y[b+i] += a[j] * x[b+ja[j]];
        }
}

void ApplyPreconditioner(void *vx, void *vy, int *blockSize, primme_params *primme)
/*
    But
    ===
    Applique le préconditionneur Multi-grid au solveur PRIMME

    Arguments
    =========
    vx          (input) - vecteur(s) d'entrée
    vy          (output) - vecteur(s) de produit vy = B^-1 * vx
    blockSize   (input) - nombre de vecteurs d'entrée
    primme      (input) - paramètres fournis par primme pour optimiser le calcul
                            (pas utilisé)

*/
{
    /* ----- DECLARATION DES VARIABLES ----- */
    // nombre de niveaux pour le multi-grid
    int totalLevels = 4;
    // nombre d'itérations pour le pre-smoothing
    int it_pre = 2;
    // nombre d'itérations pour le post-smoothing
    int it_post = 2;
    // tolérance admise pour définir le bord de la mebrane
    double tol = 1e-4;
    // niveau actuel de l'algorithme multi-grid
    int level;
    // rafinement de la grille la plus grossière
    int mCoarse = (m + pow(2, totalLevels-1) - 1) / pow(2, totalLevels-1);
    // tableaux de pointeurs pour stocker les tableaux de tous les niveaux
    int **iaM, **jaM, **born_coordM, *nM;
    double **aM, **bM, **uM;
    int i, b;
    double *x = vx, *y=vy;
    
    // on vérifie que le problème est bien posé avec les paramètres m et levels
    if(checkMultigrid(m, L, born_m, nR, tol, totalLevels)) { }
    
    /* ----- ALOCATION DES TABLEAUX ----- */
    iaM = malloc(totalLevels * sizeof(int*));
    jaM = malloc(totalLevels * sizeof(int*));
    aM = malloc(totalLevels * sizeof(double*));
    bM = malloc(totalLevels * sizeof(double*));
    uM = malloc(totalLevels * sizeof(double*));
    born_coordM = malloc(totalLevels * sizeof(int*));
    nM = malloc(totalLevels * sizeof(int));
    if (iaM == NULL || jaM == NULL || aM == NULL || bM == NULL || uM == NULL || born_coordM == NULL || nM == NULL) {
        printf("\n ERREUR : pas de mémoire pour la méthode multi-grid\n\n");
    }
    
    
    /* ----- INITIALISATION ---- */
    init_multiGrid(totalLevels, nR, L, m, born_m, &iaM, &jaM, &aM, &bM, &uM, &born_coordM, &nM);
    
    /* ----- MULTI-GRID METHOD ----- */
    for(b = 0; b < (*blockSize)*n; b+=n) {
        for(i = 0; i < n; i++) {
            bM[0][i] = x[b+i];
        }
        level = -1;
        multiGrid(it_pre, it_post, uM, iaM, jaM, aM, bM, born_coordM, nM, m, nR, L, mCoarse, level, 1);
        
        for(i = 0; i < n; i++) {
            y[b+i] = uM[0][i];
        }
    }
    
    free(iaM); free(jaM); free(aM); free(bM); free(uM);
    free(born_coordM); free(nM);
}

/**/

int primme(int primme_n, int *primme_ia, int *primme_ja, double *primme_a, 
           int nev, double *evals, double *evecs, int largest, int primme_nR, double primme_L, int primme_m, double *primme_born_m)

/*
   But
   ===
   (1) si largest = 0 :
   Calcule les nve valeurs propres les plus basses de la matrice A stoquée 
   dans le format CSR à l'aide du scalaire primme_n et des vecteurs 
   primme_ia, primme_ja et primme_ia.
 
   (2) si largest = 1 :
   calcul la plus grande valeur propres et son vecteur propre
   correspondant de la matrice A stoquée
   dans le format CSR à l'aide du scalaire primme_n et des vecteurs
   primme_ia, primme_ja et primme_ia.

  Arguments
  =========
  primme_n  (input) - le nombre d'inconus dans le système
  primme_ia (input) - le tableau 'ia' de la matrice A
  primme_ja (input) - le tableau 'ja' de la matrice A
  primme_a  (input) - le tableau 'a' de la matrice A
  nev       (input) - le nombre de valeurs propres recherchées
  evals    (output) - le tableau des valeurs propres
  evecs    (output) - le tableau des vecteurs propres
  largest  (input) - option relative au type de valeur propre calculée
                      -> 1 = la plus grande valeur propre est calculée
                      -> !1 = la plus petite valeur propre est calculée

  Retourne 0 si le caclul s'est bien déroulé, 1 si non.
*/
{
    printf("======================================================\n");
    int err;

    // norme des résidus
    double *resn = malloc(nev * sizeof(double));
    if (resn == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour un vecteur auxilière dans la fonction primme\n\n");
        return 1;
    }

    // sauvgarder les pointeurs dans des variables statiques
    n = primme_n;
    a = primme_a;
    ja = primme_ja;
    ia = primme_ia;
    nR = primme_nR;
    L = primme_L;
    m = primme_m;
    born_m = primme_born_m;

    // encoder les paramètres de PRIMME
    primme_params primme;
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; // MV product
    primme.n = primme_n;         // Matrix dimensions
    primme.numEvals = nev; // Number of wanted eigenpairs
    // primme.eps = 0.003;
    primme.printLevel = 3; // 1-4
    if(largest) {
        // calcul la plus grande valeur propre de A
        primme.target = primme_largest;
    } else {
        // calcul la plus petite valeur propre de A
        primme.target = primme_smallest;
    }
    primme.applyPreconditioner = ApplyPreconditioner;
    primme.correctionParams.precondition = 1;

    if(err = primme_set_method (DEFAULT_MIN_TIME, &primme)){
        printf("\nPRIMME: erreur N % dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
  
    // afficher les papramètres de PRIMME
    primme_display_params (primme);

    // Caclul des valeurs et vecteurs propres
    if(err = dprimme (evals, evecs, resn, &primme)){
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }

    printf("======================================================\n");
    
    // libérer la mémoire 
    primme_Free (&primme); free(resn);

    return 0;
}

