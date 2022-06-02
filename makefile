EXECS=zadanie2
MPICC?=g++

all: ${EXECS}

zadanie2: zadanie2.cpp
	${MPICC} -fopenmp zadanie2.cpp -o zadanie2 

clean:
	rm -f ${EXECS}
