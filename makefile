zadanie2:
	CC -O2 -fopenmp -std=c++2a --std=c++2a zadanie2.cpp data_source.cpp -o zadanie2

clean :
	rm -f *.o *.out *.err zadanie2
