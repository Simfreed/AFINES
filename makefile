#
## TODO: Move `libmongoclient.a` to /usr/local/lib so this can work on production servers
 
CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
BUILDDIR_DEBUG := build_debug
TARGET := bin/nt
#  
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_DEBUG := $(patsubst $(SRCDIR)/%,$(BUILDDIR_DEBUG)/%,$(SOURCES:.$(SRCEXT)=.o))

CFLAGS := -O3 -Wall -Wno-missing-braces -std=c++11 -DBOOST_TEST_DYN_LINK # -fopenmp # -g
CFLAGS_DEBUG := -Wall -std=c++11 -DBOOST_TEST_DYN_LINK -pg # -fopenmp -pg

#BOOSTSUFFIX := "-mt"
LIB := -L ${BOOST_ROOT} -lboost_unit_test_framework${BOOST_SUFFIX} -lboost_program_options${BOOST_SUFFIX}
# FOR MIDWAY
# LIB := -L /opt/local/lib/ -lboost_unit_test_framework-mt -lboost_program_options-mt # FOR MAC

INC := -I include  -I /usr/include/ -I /usr/local/include/ -I /opt/local/include/

NOW := $(shell date +"%c" | tr ' :' '_')

$(TARGET): $(OBJECTS)
	  @echo " Linking..."
	    @echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR)
	    @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR_DEBUG)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR_DEBUG)
	    @echo " $(CC) $(CFLAGS_DEBUG) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS_DEBUG) $(INC) -c -o $@ $<

clean:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

clean_debug:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR_DEBUG) $(TARGET)"; $(RM) -r $(BUILDDIR_DEBUG) $(TARGET)

tar:
	tar cfv tars/amxbd.tar src/*.cpp include/*.h

# Programs
persistence_length: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/persistence_length.cpp $(INC) $(LIB) -o bin/pl
network: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network.cpp $(INC) $(LIB) -o bin/nt
filament_force_extension: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/filament_force_extension.cpp $(INC) $(LIB) -o bin/ffe
network_pull: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network_pull.cpp $(INC) $(LIB) -o bin/ntp
debug: $(OBJECTS_DEBUG)
	$(CC) $(CFLAGS_DEBUG) $(OBJECTS_DEBUG) prog/network.cpp $(INC) $(LIB) -o bin/nt_debug
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

test:actin_tester link_tester filament_tester motor_tester # filament_ensemble_tester motor_ensemble_tester

.PHONY: clean

