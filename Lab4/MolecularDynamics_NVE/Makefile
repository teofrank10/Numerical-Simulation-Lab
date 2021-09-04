randomdir:= ../../rnumgen
seeddir:= ../../seed

CPP:=g++
CPPFLAGS:=-Wall -Wextra -pedantic
LD:=g++
CXXFLAGS := -Wall -O3 -Wextra -pedantic --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=

objdir/MolDyn_NVE.o: MolDyn_NVE.cpp MolDyn_NVE.h ${randomdir}/random.h
	mkdir -p objdir
	${CPP} ${CXXFLAGS} -c MolDyn_NVE.cpp -o objdir/MolDyn_NVE.o
