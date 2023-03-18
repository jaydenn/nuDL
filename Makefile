CC = g++
FLAGS = -O3 -Wall -g

LIBS = -L/usr/local/lib -lgsl -lgslcblas -lstdc++ -lm
INCLUDE = -I./source/include
OBJECTS = source/nuDL.o source/fileIO.o source/calcRates.o source/interactiveInput.o source/confInterval.o source/discLimit.o source/exclusionLimit.o source/nuFlux.o source/detectorFunctions.o source/formfactorSI.o source/monteCarlo.o source/likelihood.o source/nuRate.o source/SMrate.o source/BSMrate.o source/sterileOscillation.o source/sterileRate.o source/NSIrate.o

default: nuDL

nuDL: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

source/%.o: source/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm source/*.o
	-rm -f ./nuDL
