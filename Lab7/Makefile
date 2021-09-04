randomdir:= ../rnumgen
seeddir:= ../seed
MCdir:= MonteCarlo_NVT

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${MCdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
.SECONDARY:

all: ex7.1 
	cd ${randomdir} && $(MAKE) 
	
ex7.1: objdir/ex7.1.o ${randomdir}/objdir/random.o ${MCdir}/Monte_Carlo_NVT.o
	$(LD) $(CXXFLAGS) objdir/ex7.1.o ${randomdir}/objdir/random.o ${MCdir}/Monte_Carlo_NVT.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h ${MCdir}/Monte_Carlo_NVT.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/$*.cpp -o $@
	

.PHONY: all, clean, reset

clean :
	rm -f objdir/*
	rm -f ex7.1
	
reset : 
	cd ${MCdir} && ./clean.sh && cd config && rm config.final






