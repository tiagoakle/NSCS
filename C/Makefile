CC=llvm-gcc-4.2
CFLAGS=-fPIC


OBJECTDIR=./build
SRCDIR=./src
CPPSRCDIR=./src/cpp
TESTDIR=./test

TESTSOURCES=test_util.c test_linear_solver.c test_smatvec.c main_tests.c test_barriers.c test_centering_measure.c test_nscs.c
TESTOBJECTS=$(TESTSOURCES:.c=.o)

CPPSOURCES=
CPPOBJECTS=$(CPPSOURCES:.cpp=.o)

NSCSSOURCES=nscs.c spmat.c Barriers.c linear_solvers.c smatvec.c line_search.c test_util.c matlab_wrappers.c eval_cent_meas.c
NSCSOBJECTS= $(NSCSSOURCES:.c=.o)

STATICLOPENBLAS=/Users/tiago/OpenBLAS/libopenblas.a

LIBDIR=-L ../../SuiteSparse/UMFPACK/Lib\
	    -L ../../SuiteSparse/COLAMD/Lib \
	    -L ../../SuiteSparse/CHOLMOD/Lib \
	    -L ../../SuiteSparse/SuiteSparse_config \
	    -L ../../SuiteSparse/AMD/Lib \
        -L ../../SuiteSparse/CSparse/Lib \
		-L ./Lib

LIBS=-lumfpack \
	  -lcolamd \
	  -lcholmod \
	  -lsuitesparseconfig \
	  -lamd \
	  -lcsparse \
	  -lopenblas

IDIRS= -I ../../SuiteSparse/AMD/Include \
       -I ../../SuiteSparse/SuiteSparse_config \
       -I ../../SuiteSparse/UMFPACK/Include \
	   -I ../../SuiteSparse/CSparse/Include \
	   -I ./Include \
	   -I ../../ \
	   -I ./test

DYLIB=nscs.dylib
LIBOUTDIR=./Lib
DYLIBFLAGS=-dynamiclib -undefined error 


all: $(DYLIB) test

$(OBJECTDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(OBJECTDIR)/%.o: $(TESTDIR)/%.c
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(DYLIB): $(addprefix $(OBJECTDIR)/,$(NSCSOBJECTS))
	$(CC) $(DYLIBFLAGS) -g $^ -o $(LIBOUTDIR)/$@ $(LIBDIR) $(LIBS)

test: $(addprefix $(OBJECTDIR)/,$(TESTOBJECTS)) $(addprefix $(OBJECTDIR)/,$(OBJECTS)) $(addprefix $(OBJECTDIR)/,$(NSCSOBJECTS)) $(addprefix $(OBJECTDIR)/,$(CPPOBJECTS)) 
	$(CC) -undefined error $^ -g -o testprog $(LIBDIR) $(LIBS) -lcheck

clean:
	rm $(OBJECTDIR)/*.o
	rm testprog
	rm $(LIBOUTDIR)/$(DYLIB)

rtest: clean all 
	./testprog
