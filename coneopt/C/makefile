CC=gcc
CFLAGS=-fPIC
LDFLAGS=-W1,-soname,liblinear_solver.so.1
LIBS=-lumfpack -lcxsparse
SOURCES=linear_solvers.c
OBJECTS=$(SOURCES:.c=.o)

LIBRARY=liblinear_solver.so.1.0
DYLIB=libmatlab_shared.dylib
DYLIBFLAGS=-dynamiclib -undefined suppress -flat_namespace 
LIBDIR=/Users/tiago/NSCS/suite_sparse_explorations/SuiteSparse/CSparse/Lib

all: $(OBJECTS) $(DYLIB)

%.o: %.c
	$(CC) $(CFLAGS) -L$(LIBDIR) -c $< -o $@

$(DYLIB): $(OBJECTS)
	$(CC) $(DYLIBFLAGS) $< -o $(DYLIB)

clean:
	rm $(OBJECTS)
	rm $(DYLIB)
