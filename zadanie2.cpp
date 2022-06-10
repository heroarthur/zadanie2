#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "sampleSortTuple2.cpp"
#include "sampleSortTuple3.cpp"
#include "rebucketing.cpp"
// #include "reorder.cpp"
// #include "shift.cpp"

#include "generateTestData.cpp"

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

	//MPI_Tuple3 
	int blockcountTuple3[3]={1,1,1};
    MPI_Aint offsetsTuple3[3] = {offsetof(Tuple3, B), offsetof(Tuple3, B2), offsetof(Tuple3, i)};
    MPI_Datatype dataTypleTuple3[3] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    MPI_Type_create_struct(3, blockcountTuple3, offsetsTuple3, dataTypleTuple3, &MPI_Tuple3);
    MPI_Type_commit(&MPI_Tuple3);

	// MPI_Tuple2
	int blockcountTuple2[2]={charArrayLen,1};
    MPI_Aint offsetsTuple2[2] = {offsetof(Tuple2, B), offsetof(Tuple2, i)};
    MPI_Datatype dataTypleTuple2[2] = {MPI_CHAR, MPI_LONG_LONG_INT};
    MPI_Type_create_struct(2, blockcountTuple2, offsetsTuple2, dataTypleTuple2, &MPI_Tuple2);
    MPI_Type_commit(&MPI_Tuple2);

	// // ISA_Data
	// int blockcountArrData[2]={1,1};
    // MPI_Aint offsetsArrData[2] = {offsetof(ISA_Data, SA_I), offsetof(ISA_Data, B_i)};
    // MPI_Datatype dataTypeArrData[2] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    // MPI_Type_create_struct(2, blockcountArrData, offsetsArrData, dataTypeArrData, &MPI_ISA_Data);
    // MPI_Type_commit(&MPI_ISA_Data);

	// // SHIFT_data 
	// int blockcountArrData2[2]={1,1};
    // MPI_Aint offsetsArrData2[2] = {offsetof(Shift_data, i), offsetof(Shift_data, B_i_h)};
    // MPI_Datatype dataTypeArrData2[2] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    // MPI_Type_create_struct(2, blockcountArrData2, offsetsArrData2, dataTypeArrData2, &MPI_SHIFT_Data);
    // MPI_Type_commit(&MPI_SHIFT_Data);


	srand (worldRank);
    int p2 = worldSize * worldSize;

	int64 size = 20;
	vector<Tuple2>* tuple2_pointer, *tuple_sampleSorted_pointer, *tmp_pointer;
	vector<Tuple2> tuple2_Arr; tuple2_Arr.resize(size);
	vector<int64> B; B.reserve(1.2 * size);
	vector<Tuple2> tuple2_sortResult; tuple2_sortResult.reserve(1.2 * size);
	vector<Tuple2> sample; sample.resize(worldSize);
	vector<Tuple2> rootSampleRecv;
	vector<Tuple2> broadcastSample; broadcastSample.resize(worldSize - 1);
	vector<int64> pivotsPositions; pivotsPositions.resize(worldSize - 1);
	tuple2_pointer = &tuple2_Arr;
	tuple_sampleSorted_pointer = &tuple2_sortResult;

	bool allSingletones;
    
	if (worldRank == root) {
        rootSampleRecv.resize(p2);
    }



	for (int i = 0; i < size; i++) {
		fillCharArray(tuple2_Arr[i].B);
		tuple2_Arr[i].i = i;
 	}


	MPI_Barrier(MPI_COMM_WORLD);


	sample_sort_MPI_tuple2(tuple2_pointer,
						   tuple_sampleSorted_pointer,
                           &sample,
                           &rootSampleRecv,
                           &broadcastSample,
						   &pivotsPositions,
                           worldRank, 
                           worldSize);
	tmp_pointer = tuple_sampleSorted_pointer;
	tuple_sampleSorted_pointer = tuple2_pointer;
	tuple2_pointer = tmp_pointer;


	if (worldRank == 0) {
		for (int i = 0; i < tuple2_pointer->size(); i++) {
			printf("tora %s\n", tuple2_pointer->data()[i].B);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (worldRank == 1) {
		for (int i = 0; i < tuple2_pointer->size(); i++) {
			printf("tora %s\n", tuple2_pointer->data()[i].B);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (worldRank == 2) {
		for (int i = 0; i < tuple2_pointer->size(); i++) {
			printf("tora %s\n", tuple2_pointer->data()[i].B);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (worldRank == 3) {
		for (int i = 0; i < tuple2_pointer->size(); i++) {
			printf("tora %s\n", tuple2_pointer->data()[i].B);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	cout<<"rank: "<<worldRank<<" "<<tuple2_pointer->size()<<endl;


	MPI_Finalize();

	return 0;
}


