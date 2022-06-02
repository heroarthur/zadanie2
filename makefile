EXECS=zadanie2
MPICC?=g++

all: ${EXECS}

zadanie2: zadanie2.cpp
	${MPICC} -fopenmp -o zadanie2 zadanie2.cpp

clean:
	rm -f ${EXECS}
