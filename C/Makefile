CC=llvm-gcc-4.2
CFLAGS=-fPIC

TESTDIR=./test
TESTSOURCES=test_util.c test_linear_solver.c test_smatvec.c main_tests.c test_barriers.c 
TESTOBJECTS=$(TESTSOURCES:.c=.o)
TESTOBJECTDIR=./build/test

SRCDIR=./src
SOURCES=linear_solvers.c nscs_sp.c smatvec.c barriers.c ../test/test_util.c 
OBJECTS=$(SOURCES:.c=.o)
OBJECTDIR=./build

NSCSDIR=./src
NSCSSOURCES=nscs.c nscs_sp.c barriers.c linear_solvers.c smatvec.c ../test/test_util.c
NSCSOBJECTS= $(NSCSSOURCES:.c=.o)

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
	   -I ./Include \
	   -I ../../ \
	   -I ./test

DYLIB=liblinear_solvers.dylib
LIBOUTDIR=./Lib
DYLIBFLAGS=-dynamiclib -undefined error 


all: $(DYLIB) nscs

$(OBJECTDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(TESTOBJECTDIR)/%.o: $(TESTDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(DYLIB): $(addprefix $(OBJECTDIR)/,$(OBJECTS))
	$(CC) $(DYLIBFLAGS) -g $^ -o $(LIBOUTDIR)/$@ $(LIBDIR) $(LIBS)

nscs: $(addprefix $(OBJECTDIR)/,$(NSCSOBJECTS))
	$(CC) $(DYLIBFLAGS) -g $^ -o $(LIBOUTDIR)/nscs.dylib $(LIBDIR) $(LIBS)

test: $(addprefix $(TESTOBJECTDIR)/,$(TESTOBJECTS)) $(DYLIB)
	$(CC) $(addprefix $(TESTOBJECTDIR)/,$(TESTOBJECTS)) -g -o testprog $(LIBDIR) $(LIBS) -llinear_solvers -lcheck

clean:
	rm $(addprefix $(OBJECTDIR)/,$(OBJECTS))
	rm $(addprefix $(OBJECTDIR)/,$(NSCSOBJECTS))
	rm $(LIBOUTDIR)/*.dylib
	rm $(addprefix $(TESTOBJECTDIR)/,$(TESTOBJECTS)) 
	rm testprog
