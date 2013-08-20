CC=gcc
CFLAGS=-fPIC
LDFLAGS=-W1,-soname,libmatlab_shared.so.1
SOURCES=test_funcs.c matlab_shared.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=hello
LIBRARY=libmatlab_shared.so.1.0.1
DYLIB=libmatlab_shared.dylib
DYLIBFLAGS=-dynamiclib -undefined suppress -flat_namespace 

all: $(OBJECTS) $(EXECUTABLE) $(LIBRARY) $(DYLIB)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(EXECUTABLE)

$(LIBRARY): matlab_shared.o
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(LIBRARY)

%.o: %.c
	$(CC) -c $< -o $@

$(DYLIB): matlab_shared.o
	$(CC) $(DYLIBFLAGS) $< -o $(DYLIB)
clean:
	rm $(OBJECTS)
	rm $(EXECUTABLE)
	rm $(LIBRARY)
	rm $(DYLIB)
