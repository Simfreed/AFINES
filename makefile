CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
BUILDDIR_DEBUG := build_debug
TARGETDIR := bin
TARGET := bin/afines

#  

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_DEBUG := $(patsubst $(SRCDIR)/%,$(BUILDDIR_DEBUG)/%,$(SOURCES:.$(SRCEXT)=.o))

CFLAGS := -O3 -Wall -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -DBOOST_TEST_DYN_LINK -march=native# -fopenmp
CFLAGS_DEBUG := -Wall -Wunused -Wunreachable-code -std=c++11 -DBOOST_TEST_DYN_LINK -pg 

# BOOST_SUFFIX := -mt
LIB := -L ${BOOST_ROOT} -lboost_unit_test_framework${BOOST_SUFFIX} -lboost_program_options${BOOST_SUFFIX} -lboost_filesystem${BOOST_SUFFIX} -lboost_system${BOOST_SUFFIX}
INC := -I include # -isystem /usr/include/  -isystem /usr/local/include/ -isystem /opt/local/include/

#NOW := $(shell date +"%c" | tr ' :' '_')

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
	tar -cvzf tars/afines.tar.gz src/*.cpp include/*.h prog/*.cpp test/*.cpp makefile

# Programs
network: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network.cpp $(INC) $(LIB) -o bin/afines
debug: $(OBJECTS_DEBUG)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS_DEBUG) $(OBJECTS_DEBUG) prog/network.cpp $(INC) $(LIB) -o bin/afines_debug
bundles: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/bundles.cpp $(INC) $(LIB) -o bin/bun
bundles_debug: $(OBJECTS_DEBUG)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS_DEBUG) $(OBJECTS_DEBUG) prog/bundles.cpp $(INC) $(LIB) -o bin/bun_debug

# THE FOLLOWING PROGRAMS MAY OR MAY NOT EXIST; CHECK YOUR PROG FOLDER
filament_force_extension: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/filament_force_extension.cpp $(INC) $(LIB) -o bin/ffe
network_pull: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/network_pull.cpp $(INC) $(LIB) -o bin/ntp
2fil: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/2fil.cpp $(INC) $(LIB) -o bin/2f
motorwalk: $(OBJECTS)
	mkdir -p $(TARGETDIR)
	$(CC) $(CFLAGS) $(OBJECTS) prog/motors_on_struct.cpp $(INC) $(LIB) -o bin/motor_walk


# Tests
bead_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/bead_test.cpp $(INC) $(LIB) -o bin/bead_tester

spring_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/spring_test.cpp $(INC) $(LIB) -o bin/spring_tester

filament_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/filament_test.cpp $(INC) $(LIB) -o bin/filament_tester

filament_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/filament_ensemble_test.cpp $(INC) $(LIB) -o bin/filament_ensemble_tester

motor_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/motor_test.cpp $(INC) $(LIB) -o bin/motor_tester

spacer_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/spacer_test.cpp $(INC) $(LIB) -o bin/spacer_tester

motor_ensemble_tester: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test/motor_ensemble_test.cpp $(INC) $(LIB) -o bin/motor_ensemble_tester

test:bead_tester spring_tester filament_tester motor_tester

.PHONY: clean
