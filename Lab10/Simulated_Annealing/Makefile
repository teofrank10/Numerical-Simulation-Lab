randomdir:= ../../rnumgen
seeddir:= ../seed


CPP:=g++
CPPFLAGS:=-Wall -Wextra -pedantic
LD:=g++
CXXFLAGS := -Wall -O3 -Wextra -pedantic --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=
	
objdir/Simulated_Annealing.o: Simulated_Annealing.cpp Simulated_Annealing.h ${randomdir}/random.h
	mkdir -p objdir
	${CPP} ${CXXFLAGS} -c Simulated_Annealing.cpp -o objdir/Simulated_Annealing.o
	
clean :
	cd objdir && rm *.o
	



