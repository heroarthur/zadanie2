EXECS=zadanie2
MPICC?=g++

all: ${EXECS}

zadanie2: zadanie2.cpp
	${MPICC} -fopenmp --std=c++2a zadanie2.cpp -o zadanie2 

clean:
	rm -f ${EXECS}
