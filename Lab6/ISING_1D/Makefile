randomdir:= ../../rnumgen
seeddir:= ../../seed

CPP:=g++
CPPFLAGS:=-Wall -Wextra -pedantic
LD:=g++
CXXFLAGS := -Wall -O3 -Wextra -pedantic --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=

objdir/Monte_Carlo_ISING_1D.o: Monte_Carlo_ISING_1D.cpp Monte_Carlo_ISING_1D.h ${randomdir}/random.h
	mkdir -p objdir
	${CPP} ${CXXFLAGS} -c Monte_Carlo_ISING_1D.cpp -o objdir/Monte_Carlo_ISING_1D.o
	
clean :
	cd objdir && rm *.o

