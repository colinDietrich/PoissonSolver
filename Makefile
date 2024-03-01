# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
#L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a

# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ 

LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L8) $(LIBP) $(LIBEVSL) -lm -lblas -llapack

COPT = -O3

default: main

clean: 
	rm *.o 
	rm main
	rm *.txt

main: main.c prob.o time.o  umfpack_interface.o check_discretisation.o conditionBord.o convert_coord.o createShape.o dirichlet.o residu.o unknown_calculator.o write.o iterative_method.o coarse_grid.o fine_grid.o printVector.o two_grid.o multiGrid.o interface_primme.o matrixOperations.o
	cc $(COPT) $^ -o $@ $(LIB)

umfpack_interface.o: umfpack_interface.c
	cc $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

%.o: %.c
	cc $(COPT) -c $< -o $@ $(INCP)