randomdir:= ../rnumgen
seeddir:= ../seed
moleculedir:= MolecularDynamics_NVE

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -I ${moleculedir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
ex4.2: objdir/ex4.2.o ${randomdir}/objdir/random.o ${moleculedir}/MolDyn_NVE.o
	$(LD) $(CXXFLAGS) objdir/ex4.2.o ${randomdir}/objdir/random.o ${moleculedir}/MolDyn_NVE.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h ${moleculedir}/MolDyn_NVE.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/$*.cpp -o $@
	
ex4.1: objdir/ex4.1.o ${randomdir}/objdir/random.o ${moleculedir}/MolDyn_NVE.o
	$(LD) $(CXXFLAGS) objdir/ex4.1.o ${randomdir}/objdir/random.o ${moleculedir}/MolDyn_NVE.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h ${moleculedir}/MolDyn_NVE.h
	mkdir -p objdir
	$(LD) $(CXXFLAGS) -c source/$*.cpp -o $@


.SECONDARY:

all: ex4.1 ex4.2
	cd ${randomdir} && $(MAKE) 
	
ex1.%: objdir/ex1.%.o ${randomdir}/objdir/random.o
	$(CC) $(CXXFLAGS) objdir/ex1.$*.o ${randomdir}/objdir/random.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h
	mkdir -p objdir
	$(CC) $(CXXFLAGS) -c source/$*.cpp -o $@

.PHONY: all, clean, reset

clean :
	rm -f objdir/*
	rm -f ex4.1 ex4.2 
	
reset : 
	cd ${moleculedir} && ./clean.sh && cd config && rm old.0 old.final






