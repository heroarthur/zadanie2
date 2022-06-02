#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "sampleSort.cpp"

#define int64 long long int

using namespace std;





// void sample_sort_MPI()




int main(int argc, char** argv) {

	// MPI_Init(&argc, &argv);

	srand (time(NULL));

	int size = 100000000;
	int64* A = (int64*) malloc (size * sizeof(int64));
	int64* A_test = (int64*) malloc (size * sizeof(int64));
	
	for (int i = 0; i < size; i++) {
		A[i] = (rand() % 10000) + 1;
		A_test[i] = A[i];
		// cout<<A[i]<<" ";
	}

	// cout<<endl<<endl;


	// local_sort_openMP(A, size);

	// std::sort(A_test, A_test + size, less<int64>());
	
	// for (int i = 0; i < size; i++) {
	// 	cout<<A[i]<<" ";
	// }
	// cout <<endl<<endl;

	// MPI_Finalize();

	return 0;
}



// real    0m37,157s
// user    0m36,686s
// sys     0m0,460s