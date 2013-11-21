CC=llvm-g++-4.2
CFLAGS=-fPIC

OBJECTDIR=./build
SRCDIR=./src
TESTDIR=./test

NSCSSOURCES= nscs_long_step.cpp
NSCSOBJECTS= $(NSCSSOURCES:.cpp=.o)

STATICLOPENBLAS=/Users/tiago/OpenBLAS/libopenblas.a

LIBDIR= -L ./Lib

LIBS= -lopenblas

IDIRS= -I ./Include 

DYLIB=nscs.dylib
LIBOUTDIR=./build
DYLIBFLAGS=-dynamiclib -undefined error 

all: $(DYLIB)

$(OBJECTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(OBJECTDIR)/%.o: $(TESTDIR)/%.cpp
	$(CC) $(CFLAGS) -g $(IDIRS) -c $< -o $@

$(DYLIB): $(addprefix $(OBJECTDIR)/,$(NSCSOBJECTS))
	$(CC) $(DYLIBFLAGS) -g $^ -o $(LIBOUTDIR)/$@ $(LIBDIR) $(LIBS)

clean:
	rm $(OBJECTDIR)/*.o
	rm $(LIBOUTDIR)/$(DYLIB)
