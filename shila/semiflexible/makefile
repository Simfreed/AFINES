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
CFLAGS := -g -Wall -std=c++11
LIB := -L lib
INC := -I include

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

test:actin_tester actin_ensemble_tester link_tester link_ensemble_tester motor_tester motor_ensemble_tester

# Spikes
ticket:
	  $(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean

