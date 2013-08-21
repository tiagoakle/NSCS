CC=llvm-gcc-4.2
CFLAGS=-fPIC
LDFLAGS=-W1,-soname,liblinear_solver.so.1

TESTDIR=./test
TESTSOURCES=test_util.c test_linear_solver.c 
TESTOBJECTS=./build/test/test_util.o ./build/test/test_linear_solver.o
TESTOBJECTDIR=./build/test

SRCDIR=./src
SOURCES=linear_solvers.c
OBJECTS=$(SOURCES:.c=.o)
OBJECTDIR=./build

STATICLOPENBLAS=/Users/tiago/OpenBLAS/libopenblas.a

LIBDIR=-L ./SuiteSparse/UMFPACK/Lib\
	    -L ./SuiteSparse/COLAMD/Lib \
	    -L ./SuiteSparse/CHOLMOD/Lib \
	    -L ./SuiteSparse/SuiteSparse_config \
	    -L ./SuiteSparse/AMD/Lib \
        -L ./SuiteSparse/CSparse/Lib \
		-L ./Lib

LIBS=-lumfpack \
	  -lcolamd \
	  -lcholmod \
	  -lsuitesparseconfig \
	  -lamd \
	  -lcsparse \
	  -lopenblas 

IDIRS= -I ./SuiteSparse/AMD/Include \
       -I ./SuiteSparse/SuiteSparse_config \
       -I ./SuiteSparse/UMFPACK/Include \
	   -I ./SuiteSparse/CSparse/Include \
	   -I ./Include

DYLIB=liblinear_solvers.dylib
LIBOUTDIR=./Lib
DYLIBFLAGS=-dynamiclib -undefined error 


all: $(DYLIB) test

$(OBJECTDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(TESTOBJECTDIR)/%.o: $(TESTDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(DYLIB): $(OBJECTDIR)/$(OBJECTS)
	$(CC) $(DYLIBFLAGS) -g $(OBJECTDIR)/*.o -o $(LIBOUTDIR)/$(DYLIB) $(LIBDIR) $(LIBS)

test: $(TESTOBJECTS) $(DYLIB)
	$(CC) $(TESTOBJECTDIR)/*.o -o testprog $(LIBDIR) $(LIBS) -llinear_solvers

clean:
	rm $(OBJECTDIR)/$(OBJECTS)
	rm $(LIBOUTDIR)/$(DYLIB)
	rm $(TESTOBJECTDIR)/$(TESTOBJECTS)
	rm testprog
