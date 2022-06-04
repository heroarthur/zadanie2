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

	int worldRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	int worldSize;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int blockcount[3]={1,1,1};
    MPI_Aint offsets[3] = {offsetof(Tuple3, B), offsetof(Tuple3, B2), offsetof(Tuple3, i)};
    MPI_Datatype dataType[3] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    
    MPI_Type_create_struct(3, blockcount, offsets, dataType, &MPI_Tuple3);
    MPI_Type_commit(&MPI_Tuple3);

	srand (worldRank);
    int p2 = worldSize * worldSize;

	int64 size = 10000000;
	vector<Tuple3>* tuple3_pointer, *tuple_sampleSorted_pointer, *tmp_pointer;
	vector<Tuple3> tuple3_Arr; tuple3_Arr.resize(size);
	vector<int64> B; B.reserve(1.5 * size);
	vector<Tuple3> tuple3_sortResult; tuple3_sortResult.reserve(1.5 * size);
	vector<Tuple3> sample; sample.resize(worldSize);
	vector<Tuple3> rootSampleRecv;
	vector<Tuple3> broadcastSample; broadcastSample.resize(worldSize - 1);
	vector<int64> pivotsPositions; pivotsPositions.resize(worldSize - 1);
	tuple3_pointer = &tuple3_Arr;
	tuple_sampleSorted_pointer = &tuple3_sortResult;

	bool allSingletones;
    
	if (worldRank == root) {
        rootSampleRecv.resize(p2);
    }



	
	for (int i = 0; i < size; i++) {
		tuple3_Arr[i].B = (rand() % 50) + 1;
		tuple3_Arr[i].B2 = (rand() % 50) + 1;
		tuple3_Arr[i].i = i;
 	}

	sample_sort_MPI_tuple3(tuple3_pointer,
						   tuple_sampleSorted_pointer,
                           &sample,
                           &rootSampleRecv,
                           &broadcastSample,
						   &pivotsPositions,
                           &size, 
                           worldRank, 
                           worldSize);
	tmp_pointer = tuple_sampleSorted_pointer;
	tuple_sampleSorted_pointer = tuple3_pointer;
	tuple3_pointer = tmp_pointer;


	B.clear();
	rebucketing_2h_group_rank(tuple3_pointer, 
                              &B, 
                              &allSingletones,
                              worldRank,
                              worldSize);



	MPI_Finalize();

	return 0;
}



// real    0m37,157s
// user    0m36,686s
// sys     0m0,460s