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
#include "reorder.cpp"
#include "shift.cpp"

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

	// TwoInts64
	int blockcountArrData[2]={1,1};
    MPI_Aint offsetsArrData[2] = {offsetof(TwoInts64, i1), offsetof(TwoInts64, i2)};
    MPI_Datatype dataTypeArrData[2] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
    MPI_Type_create_struct(2, blockcountArrData, offsetsArrData, dataTypeArrData, &MPI_TwoInts64);
    MPI_Type_commit(&MPI_TwoInts64);



	srand (worldRank);
    int p2 = worldSize * worldSize;

	int64 singleNodeDataSize = 100;
	vector<Tuple2> *tuple2_pointer, *tuple2_help_pointer, *tmp_pointer;
	vector<Tuple2> tuple2_Arr; tuple2_Arr.reserve(1.2 * singleNodeDataSize); tuple2_Arr.resize(singleNodeDataSize);
	
	vector<int64> B_1; B_1.reserve(1.2 * singleNodeDataSize);
	vector<int64> B_2; B_2.reserve(1.2 * singleNodeDataSize);
	vector<int64> *B_pointer, *B_help_pointer, *B_tmp_pointer;

	vector<Tuple2> tuple2_sortResult; tuple2_sortResult.reserve(1.2 * singleNodeDataSize);
	
	HelpingVectorsSendingOperations helpVectorsSendingOperations;
	HelpingVectorsSampleSort2 helpVectorsSampleSort2;

	initializeHelpingVectorsSendingOperations(&helpVectorsSendingOperations, singleNodeDataSize, worldSize);
	initializeHelpingVectorsSampleSort2(&helpVectorsSampleSort2, worldSize);


	tuple2_pointer = &tuple2_Arr;
	tuple2_help_pointer = &tuple2_sortResult;

	B_pointer = &B_1;
	B_help_pointer = &B_2;

	vector<int64> SA;

	bool allSingletones;
    

	for (int i = 0; i < singleNodeDataSize; i++) {
		fillCharArray(tuple2_Arr[i].B);
		tuple2_Arr[i].i = i + worldRank * singleNodeDataSize;
 	}

	if (worldRank == root) {
		tuple2_Arr[0].B[0] = '$';
	}


	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < 100; i++) {
		sample_sort_MPI_tuple2(tuple2_pointer,
						tuple2_help_pointer,
						&helpVectorsSampleSort2,
						worldRank, 
						worldSize);
		switchPointersTuple2(&tuple2_pointer, &tuple2_help_pointer);

		// if (worldRank == root)
		// 	cout<<tuple2_pointer->size()<<endl;
	}




	// assign_h_group_rank(tuple2_pointer, 
	// 				B_pointer, 
	// 				worldRank,
	// 				worldSize);
	
	// print_MPI_vector(B_pointer, worldRank, worldSize);

	// initialize_SA(&SA, tuple2_pointer);

	// for (int i = 0; i < 50; i++) {
	// 	reorder_and_rebalance(&B_pointer, 
	// 	                      &B_help_pointer, 
	// 	                      &SA,
	// 						  &helpVectorsSendingOperations,
	// 	                      worldRank, 
	// 	                      worldSize);

	// 	shift_by_h(&B_pointer, 
	// 			   &B_help_pointer, 
	// 			   &SA,
	// 			   &helpVectorsSendingOperations,
	// 			   10,
	// 			   worldRank, 
	// 			   worldSize);
	// }



	// reorder_and_rebalance(&B_pointer, 
    //                       &B_help_pointer, 
    //                       &SA,
	// 					  	 &dataForPartitions,
    //                       worldRank, 
    //                       worldSize);



	// // cout<<"rank "<<worldRank<<" "<<B_pointer->size()<<endl;
	// // print_MPI_vector(B_pointer, worldRank, worldSize);

	// for (int i = 0; i < 5000; i++) {
	// 	shift_by_h(&B_pointer, 
	// 			&B_help_pointer, 
	// 			&SA,
	// 			&dataForPartitions,
	// 			10,
	// 			worldRank, 
	// 			worldSize);
	// }



	// print_MPI_vector(B_pointer, worldRank, worldSize);

	MPI_Finalize();

	return 0;
}


