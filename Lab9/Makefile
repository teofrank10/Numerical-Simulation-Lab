randomdir:= ../rnumgen
seeddir:= ../seed
gadir:= Genetic_Algorithm

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${gadir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
.SECONDARY:

all: ex9.1 
	cd ${randomdir} && $(MAKE) 
	
ex9.1: objdir/ex9.1.o ${randomdir}/objdir/random.o ${gadir}/objdir/Genetic_Algorithm.o
	$(LD) $(CXXFLAGS) objdir/ex9.1.o ${randomdir}/objdir/random.o ${gadir}/objdir/Genetic_Algorithm.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h ${gadir}/Genetic_Algorithm.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/$*.cpp -o $@

.PHONY: all, clean, reset

clean :
	rm -f objdir/*
	rm -f ex9.1

	



