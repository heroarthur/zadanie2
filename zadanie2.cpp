#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "sampleSort.cpp"
#include "rebucketing.cpp"

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

	srand (wordRank);
    int p2 = worldSize * worldSize;

	int64 size = 1000;
	vector<Tuple3>* A_pointer, *A_sampleSorted_pointer, *tmp_pointer;
	vector<Tuple3> A; A.resize(size);
	vector<Tuple3> A_sampleSorted; A_sampleSorted.reserve(1.5 * size);
	vector<Tuple3> sample; sample.resize(worldSize);
	vector<Tuple3> rootSampleRecv;
	vector<Tuple3> broadcastSample; broadcastSample.resize(worldSize - 1);
	vector<int64> pivotsPositions; pivotsPositions.resize(worldSize - 1);
	A_pointer = &A;
	A_sampleSorted_pointer = &A_sampleSorted;

    
	if (wordRank == root) {
        rootSampleRecv.resize(p2);
    }



	
	for (int i = 0; i < size; i++) {
		A[i].B = (rand() % 50) + 1;
		A[i].B2 = (rand() % 50) + 1;
	}

	sample_sort_MPI_tuple3(A_pointer,
						   A_sampleSorted_pointer,
                           &sample,
                           &rootSampleRecv,
                           &broadcastSample,
						   &pivotsPositions,
                           &size, 
                           wordRank, 
                           worldSize);
	tmp_pointer = A_sampleSorted_pointer;
	A_sampleSorted_pointer = A_pointer;
	A_pointer = tmp_pointer;


	MPI_Finalize();

	return 0;
}



// real    0m37,157s
// user    0m36,686s
// sys     0m0,460s