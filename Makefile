CC=g++
CFLAGS=-O2 -march=x86-64 -ffast-math -fopenmp
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L/Your_library
INCLUDE=-I/Your_library
#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm

OBJECTS=BicMix.o

main: $(OBJECTS)
	g++ $(CFLAGS) $(OBJECTS) -o BicMix $(LINKER) $(LIBS)
BicMix.o: BicMix.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c BicMix.cpp
clean:
	rm *.o
