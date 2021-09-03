randomdir:= ../rnumgen
seeddir:= ../seed
sadir:= Simulated_Annealing
gadir:= Genetic_Algorithm

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${sadir} -D 'SEEDDIR="${seeddir}"'
CPPFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${gadir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
MPI_CPP = mpicxx
MPI_LD = mpicxx
	
.SECONDARY:

all: ex10.1 ex10.2
	cd ${randomdir} && $(MAKE) 
	
ex10.1: objdir/ex10.1.o ${randomdir}/objdir/random.o ${sadir}/objdir/Simulated_Annealing.o
	$(LD) $(CXXFLAGS) objdir/ex10.1.o ${randomdir}/objdir/random.o ${sadir}/objdir/Simulated_Annealing.o -o $@
	
ex10.2: objdir/ex10.2.o ${randomdir}/objdir/random.o ${gadir}/objdir/Genetic_Algorithm.o
	$(MPI_LD) $(CPPFLAGS) objdir/ex10.2.o ${randomdir}/objdir/random.o ${gadir}/objdir/Genetic_Algorithm.o -o $@
	
objdir/ex10.1.o: source/ex10.1.cpp ${randomdir}/random.h ${sadir}/Simulated_Annealing.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/ex10.1.cpp -o $@
	
objdir/ex10.2.o: source/ex10.2.cpp ${randomdir}/random.h ${gadir}/Genetic_Algorithm.h
	mkdir -p objdir
	$(MPI_CPP) $(CPPFLAGS) -c source/ex10.2.cpp -o $@

.PHONY: all, clean, reset

clean :
	rm -f objdir/*
	rm -f ex10.1
	rm -f ex10.2

	



