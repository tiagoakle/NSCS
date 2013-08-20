CC=llvm-gcc-4.2
CFLAGS=-fPIC
LDFLAGS=-W1,-soname,liblinear_solver.so.1
OBJECTDIR=./build
TESTOBJECTDIR=./build/test
SRCDIR=./src
SOURCES=linear_solvers.c
TESTSOURCES=test_linear_solver.c
TESTDIR=./test
OBJECTS=$(SOURCES:.c=.o)
TESTOBJECTS=$(TESTSOURCES:.c=.o)

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

test: $(TESTOBJECTDIR)/$(TESTOBJECTS) $(DYLIB)
	$(CC) $(TESTOBJECTDIR)/*.o -o testprog $(LIBDIR) $(LIBS) -llinear_solvers

clean:
	rm $(OBJECTDIR)/$(OBJECTS)
	rm $(LIBOUTDIR)/$(DYLIB)
	rm $(TESTOBJECTDIR)/$(TESTOBJECTS)
	rm testprog
