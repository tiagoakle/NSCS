LIBDIR=./SuiteSparse/CSparse/Lib 
ULIBDIR=-L ./SuiteSparse/UMFPACK/Lib\
	    -L ./SuiteSparse/COLAMD/Lib \
	    -L ./SuiteSparse/CHOLMOD/Lib \
	    -L ./SuiteSparse/SuiteSparse_config \
	    -L ./SuiteSparse/AMD/Lib

ULIBS=-lumfpack \
	  -lcolamd \
	  -lcholmod \
	  -lsuitesparseconfig \
	  -lamd \
	  -lblas

UIDIRS=-I ./SuiteSparse/AMD/Include \
       -I ./SuiteSparse/SuiteSparse_config
LIBS=-lcsparse 
EXECUTABLE=test_prog

all: test umftest

test:
	gcc test_triplet_to_ccf.c -L$(LIBDIR) $(LIBS) -o $(EXECUTABLE)

umftest:
	gcc test_umfpack.c $(UIDIRS) $(ULIBDIR) $(ULIBS) -o test_umfpack

clean:
	rm $(EXECUTABLE)
	rm test_umfpack
