CC = g++
FLAGS = #-g -O3 #-DMPI -DOMPI_SKIP_MPICXX
NESTLIBDIR = ../MultiNest_v3.9

LIBS = -L$(NESTLIBDIR) -lnest3 -llapack -lgsl -lgslcblas -lstdc++ -lm #-lgfortran
INCLUDE = -I./source/include
OBJECTS = source/Conudl.o source/nuFlux.o source/detectorFunctions.o source/detectors.o source/formfactorSI.o source/monteCarlo.o source/likelihood.o source/CNNSrate.o

default: Conudl

Conudl: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

source/%.o: source/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm source/*.o
	-rm -f ./BSMNu
