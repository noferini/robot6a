BIN=./bin

# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler, linker, ... flags:
CC_FLAGS := $(shell root-config --cflags)
LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --glibs) -lGeom -lEG -lRGL

INCLUDES := -I$(shell root-config --incdir)


START := $(shell rm -rf LinkDef.* LinkDef2.*)

CPP_FILES := $(wildcard *.cpp)
HEADERS :=  $(wildcard EEE*.h)


SRC_MAIN=prog.C
OBJ_FILES := $(addprefix ./,$(notdir $(CPP_FILES:.cpp=.o)))
OBJ_MAIN=$(SRC_MAIN:.C=.o)

# the build target executable:
MAIN = $(BIN)/prog

# build the executable (remove LinkDef at the end to avoid conflict in the next "make" call)
all:    main

# build reco executable
main:   $(OBJ_FILES) $(OBJ_MAIN) prog.C
	$(CC) $(CC_FLAGS) $(INCLUDES) -o $(MAIN) $(OBJ_FILES) prog.C $(LFLAGS) $(LIBS)

#compilation commands
.cpp.o:
	$(CC) $(CC_FLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<

.C.o:
	$(CC) $(CC_FLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<


# run the first time to build the structure
config:
	@rm -rf bin
	@echo
	@echo "do \"make\""
	@echo
	@echo "set your installation dir \"export PATH_INSTALL=...\""
	@echo
	@echo "do \"make install\" (or, if needed, \"sudo make install\")"
	@echo
	@mkdir bin
	@bash checklib

#clean commands

clean:
	rm -rf LinkDef.* *.o bin/*.exe *~

install:
	cp $(BIN)/*.exe $(PATH_INSTALL)
