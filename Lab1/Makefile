randomdir:= ../rnumgen
seeddir:= ../seed

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
.SECONDARY:

all: ex1.1 ex1.2 ex1.3
	cd ${randomdir} && $(MAKE) 
	
ex1.%: objdir/ex1.%.o ${randomdir}/objdir/random.o
	$(CC) $(CXXFLAGS) objdir/ex1.$*.o ${randomdir}/objdir/random.o -o $@
	
objdir/%.o: source/%.cpp ${randomdir}/random.h
	mkdir -p objdir
	$(CC) $(CXXFLAGS) -c source/$*.cpp -o $@

.PHONY: all, clean 

clean :
	rm -f objdir/*
	rm -f ex1.1 ex1.2 ex1.3
