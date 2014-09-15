CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=globals.cpp actin.cpp actin_ensemble.cpp Link.cpp link_ensemble.cpp motor.cpp motor_ensemble.cpp persistence_length.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=per

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


