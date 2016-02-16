#
## TODO: Move `libmongoclient.a` to /usr/local/lib so this can work on production servers
 
CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
TARGET := bin/runner
#  
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
# CFLAGS := -O3 -Wall -std=c++0x -DBOOST_TEST_DYN_LINK -openmp # -pg
CFLAGS := -Wall -std=c++11 -DBOOST_TEST_DYN_LINK -fopenmp -g
LIB := -L /software/boost-1.50-el6-x86_64/lib/ -L lib -lboost_unit_test_framework -lboost_program_options
INC := -I include  -I /usr/include/ -I /usr/local/include/

$(TARGET): $(OBJECTS)
	  @echo " Linking..."
	    @echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR)
	    @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Programs
persistence_length: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/persistence_length.cpp $(INC) $(LIB) -o bin/pl
network: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network.cpp $(INC) $(LIB) -o bin/nt
filament_force_extension: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/filament_force_extension.cpp $(INC) $(LIB) -o bin/ffe
network_pull: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network_pull.cpp $(INC) $(LIB) -o bin/ntp
debug: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network.cpp $(INC) $(LIB) -o bin/nt_debug
baoab_network: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/baoab_network.cpp $(INC) $(LIB) -o bin/baoab_nt
llf_network: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/llf_network.cpp $(INC) $(LIB) -o bin/llf_nt
2fil: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/2fil.cpp $(INC) $(LIB) -o bin/2f

# Tests
actin_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/actin_test.cpp $(INC) $(LIB) -o bin/actin_tester

actin_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/actin_ensemble_test.cpp $(INC) $(LIB) -o bin/actin_ensemble_tester

link_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/link_test.cpp $(INC) $(LIB) -o bin/link_tester

link_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/link_ensemble_test.cpp $(INC) $(LIB) -o bin/link_ensemble_tester

motor_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/motor_test.cpp $(INC) $(LIB) -o bin/motor_tester

motor_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/motor_ensemble_test.cpp $(INC) $(LIB) -o bin/motor_ensemble_tester

filament_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/filament_test.cpp $(INC) $(LIB) -o bin/filament_tester

lammps_filament_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/lammps_filament_test.cpp $(INC) $(LIB) -o bin/lammps_filament_tester

DLfilament_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/DLfilament_test.cpp $(INC) $(LIB) -o bin/DLfilament_tester

filament_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/filament_ensemble_test.cpp $(INC) $(LIB) -o bin/filament_ensemble_tester

test:actin_tester actin_ensemble_tester link_tester link_ensemble_tester motor_tester motor_ensemble_tester

# Spikes
ticket:
	  $(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean

