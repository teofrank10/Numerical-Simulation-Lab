randomdir:= ../rnumgen
seeddir:= ../seed

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=

	
.SECONDARY:

all: ex2.1 ex2.2
	cd ${randomdir} && $(MAKE) 
	
ex2.%: objdir/ex2.%.o ${randomdir}/objdir/random.o
	$(CC) $(CXXFLAGS) objdir/ex2.$*.o ${randomdir}/objdir/random.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h
	mkdir -p objdir
	$(CC) $(CXXFLAGS) -c source/$*.cpp -o $@

.PHONY: all, clean 

clean :
	rm -f objdir/*
	rm -f ex2.1 ex2.2
