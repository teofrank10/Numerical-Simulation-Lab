randomdir:= ../rnumgen
seeddir:= ../seed
isingdir:= ISING_1D

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${isingdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
.SECONDARY:

all: ex6.1 
	cd ${randomdir} && $(MAKE) 
	
ex6.1: objdir/ex6.1.o ${randomdir}/objdir/random.o ${isingdir}/Monte_Carlo_ISING_1D.o
	$(LD) $(CXXFLAGS) objdir/ex6.1.o ${randomdir}/objdir/random.o ${isingdir}/Monte_Carlo_ISING_1D.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h ${isingdir}/Monte_Carlo_ISING_1D.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/$*.cpp -o $@
	

.PHONY: all, clean, reset

clean :
	rm -f objdir/*
	rm -f ex6.1
	
reset : 
	cd ${isingdir} && ./clean.sh && cd config && rm old.0 old.final






