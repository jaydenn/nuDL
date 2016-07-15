CC = g++
FLAGS = -g #-Ofast

LIBS = -lgsl -lgslcblas -lstdc++ -lm
INCLUDE = -I./source/include
OBJECTS = source/Conudl.o source/fileIO.o source/calcRates.o source/interactiveInput.o source/discLimit.o source/nuFlux.o source/detectorFunctions.o source/formfactorSI.o source/monteCarlo.o source/likelihood.o source/nuRate.o source/SMrate.o source/BSMrate.o

default: Conudl

Conudl: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

source/%.o: source/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm source/*.o
	-rm -f ./Conudl
