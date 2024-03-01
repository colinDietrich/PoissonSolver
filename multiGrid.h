int multiGrid(int it_pre, int it_post, double **uM, int **iaM, int **jaM, double **aM, double **bM, int **born_coordM, int *nM, int m, int nR, double L, double mCoarse, int niv, int iter);

int init_multiGrid(int levels, int nR, int L, int m, double *born_m, int ***iaM, int ***jaM, double ***aM, double ***bM, double ***uM, int ***born_coordM, int **nM);
