randomdir:= ../../rnumgen
seeddir:= ../../seed

CPP:=g++
CPPFLAGS:=-Wall -Wextra -pedantic
LD:=g++
CXXFLAGS := -Wall -O3 -Wextra -pedantic --std=c++11 -I ${randomdir} -D 'SEEDDIR="${seeddir}"'
LDFLAGS:=

objdir/Monte_Carlo_NVT.o: Monte_Carlo_NVT.cpp Monte_Carlo_NVT.h ${randomdir}/random.h
	mkdir -p objdir
	${CPP} ${CXXFLAGS} -c Monte_Carlo_NVT.cpp -o objdir/Monte_Carlo_NVT.o
