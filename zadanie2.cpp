#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "sampleSort.cpp"

#define int64 long long int
#define root 0

using namespace std;





// void sample_sort_MPI()

// struct Tuple2 {
//     char B[K];
//     int64 i;
// };

// struct Tuple3 {
//     int64 B;
//     int64 B2;
//     int64 i;
// };


int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	int wordRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &wordRank);
	int worldSize;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int blockcount[3]={1,1,1};
    MPI_Aint offsets[3] = {offsetof(Tuple3, B), offsetof(Tuple3, B2), offsetof(Tuple3, i)};
    MPI_Datatype dataType[3] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    
    MPI_Type_create_struct(3, blockcount, offsets, dataType, &MPI_Tuple3);
    MPI_Type_commit(&MPI_Tuple3);

	srand (time( NULL ) * wordRank);
    int p2 = worldSize * worldSize;

	int64 size = 100;
	Tuple3* A = (Tuple3*) malloc (size * sizeof(Tuple3));
	Tuple3* sample = (Tuple3*) malloc (worldSize * sizeof(Tuple3));
	Tuple3* rootSampleRecv;
	Tuple3* broadcastSample = (Tuple3*) malloc ((worldSize - 1) * sizeof(Tuple3));
    
	if (wordRank == root) {
        rootSampleRecv = (Tuple3*) malloc (p2 * sizeof(Tuple3));
		memset (rootSampleRecv,0,p2 * sizeof(Tuple3));
    }



	
	for (int i = 0; i < size; i++) {
		A[i].B = (rand() % 50) + 1;
		A[i].B2 = (rand() % 50) + 1;
	}

	sample_sort_MPI_tuple3(A, 
                           sample,
                           rootSampleRecv,
                           broadcastSample,
                           &size, 
                           wordRank, 
                           worldSize);

	// cout<<endl<<endl;


	// local_sort_openMP(A, size);

	// std::sort(A_test, A_test + size, less<int64>());
	
	// for (int i = 0; i < size; i++) {
	// 	cout<<A[i]<<" ";
	// }
	// cout <<endl<<endl;

	MPI_Finalize();

	return 0;
}



// real    0m37,157s
// user    0m36,686s
// sys     0m0,460s