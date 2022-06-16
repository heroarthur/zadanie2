#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "SA.cpp"

using namespace std;




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

	// TwoInts64
	int blockcountArrData[2]={1,1};
    MPI_Aint offsetsArrData[2] = {offsetof(TwoInts64, i1), offsetof(TwoInts64, i2)};
    MPI_Datatype dataTypeArrData[2] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    MPI_Type_create_struct(2, blockcountArrData, offsetsArrData, dataTypeArrData, &MPI_TwoInts64);
    MPI_Type_commit(&MPI_TwoInts64);


	srand (worldRank);

    int p2 = worldSize * worldSize;
	int64 allDataSize = 100;
	int64 singleNodeDataSize = allDataSize / worldSize;
	
	vector<int64> B_1; //B_1.reserve(1.2 * singleNodeDataSize);
	vector<int64> B_2; //B_2.reserve(1.2 * singleNodeDataSize);
	// vector<int64> *B_pointer, *B_second_pointer, *B_ISA_pointer;
	vector<int64> *B_ISA_pointer, *SA_pointer;


	// vector<Tuple2> *tuple2_pointer, *tuple2_second_pointer;
	vector<Tuple2> tuple2_Arr; tuple2_Arr.resize(singleNodeDataSize); //tuple2_Arr.reserve(1.2 * singleNodeDataSize);
	vector<Tuple2> tuple2_second; //tuple2_second.reserve(1.2 * singleNodeDataSize);
	// tuple2_pointer = &tuple2_Arr;
	// tuple2_second_pointer = &tuple2_second;


	// vector<Tuple3> *tuple3_pointer, *tuple3_second_pointer;
	vector<Tuple3> tuple3; //tuple3.reserve(1.2 * singleNodeDataSize);
	vector<Tuple3> tuple3_second; //tuple3_second.reserve(1.2 * singleNodeDataSize);
	// tuple3_pointer = &tuple3;
	// tuple3_second_pointer = &tuple3_second;

	HelpingVectorsSendingOperations helpVectorsSendingOperations;
	HelpingVectorsSampleSort2 helpVectorsSampleSort2;
	HelpingVectorsSampleSort3 helpVectorsSampleSort3;

	initializeHelpingVectorsSendingOperations(&helpVectorsSendingOperations, singleNodeDataSize, worldSize);
	initializeHelpingVectorsSampleSort2(&helpVectorsSampleSort2, worldSize);
	initializeHelpingVectorsSampleSort3(&helpVectorsSampleSort3, worldSize);


	// B_pointer = &B_1;
	// B_second_pointer = &B_2;
	// B_ISA_pointer = &B_ISA;

	vector<int64> SA, SA_second;
	// vector<int64> *SA_pointer, *SA_second_pointer;
	// SA_pointer = &SA;
	// SA_second_pointer = &SA_second;

	// vector<int64> *ISA;

	bool allSingletones;

	for (int i = 0; i < 100; i++) {
		for (int i = 0; i < singleNodeDataSize; i++) {
			fillCharArray(tuple2_Arr[i].B);
			tuple2_Arr[i].i = i + worldRank * singleNodeDataSize;
		}

		if (worldRank == root) {
			tuple2_Arr[0].B[0] = '$';
		}

		
		SA_algorithm(&B_1,
					&B_2,
					&SA,
					&SA_second,
					&tuple2_Arr, 
					&tuple2_second, 
					&tuple3, 
					&tuple3_second,
					&helpVectorsSendingOperations,
					&helpVectorsSampleSort2,
					&helpVectorsSampleSort3,
					&B_ISA_pointer, 
					worldRank,
					worldSize);
		
		// print_MPI_vector(&SA, worldRank, worldSize);
	}




	MPI_Finalize();

	return 0;
}


