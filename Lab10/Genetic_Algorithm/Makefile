randomdir:= ../../rnumgen
seeddir:= ../seed


CPP:=g++
CPPFLAGS:=-Wall -Wextra -pedantic
LD:=g++
CXXFLAGS := -Wall -O3 -Wextra -pedantic --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
MPI_CPP = mpicxx
MPI_LD = mpicxx
	
objdir/Genetic_Algorithm.o: Genetic_Algorithm.cpp Genetic_Algorithm.h ${randomdir}/random.h
	mkdir -p objdir
	${MPI_CPP} ${CXXFLAGS} -c Genetic_Algorithm.cpp -o objdir/Genetic_Algorithm.o
	
clean :
	cd objdir && rm *.o
	



