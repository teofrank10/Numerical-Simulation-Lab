randomdir:= ../rnumgen
seeddir:= ../seed

CC := g++
LD:=g++
CXXFLAGS := -Wall -O3 --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=

	
.SECONDARY:

all: ex8.2 ex8.1
	cd ${randomdir} && $(MAKE) 
	
ex8.%: objdir/ex8.%.o ${randomdir}/objdir/random.o
	$(CC) $(CXXFLAGS) objdir/ex8.$*.o ${randomdir}/objdir/random.o -o $@
	
	
objdir/%.o: source/%.cpp ${randomdir}/random.h
	mkdir -p objdir
	$(CC) $(CXXFLAGS) -c source/$*.cpp -o $@

.PHONY: all, clean 

clean :
	rm -f objdir/*
	rm -f ex8.2 ex8.1
