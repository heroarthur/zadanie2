#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <numeric>


#include <omp.h>
#include <mpi.h>

#ifndef   sortsOpenMP
#define   sortsOpenMP

#include "sortsOpenMP.cpp"

#endif

#include<bits/stdc++.h>

#define root 0

#define wyslijRaz (2147483647 / 10)

using namespace std;

// struct Tuple2 {
//     char B[K];
//     int64 i;
// };

// struct Tuple3 {
//     int64 B;
//     int64 B2;
//     int64 i;
// };



int64 binarySearchTuple3(vector<Tuple3>* arr, Tuple3 tuple, int64 l, int64 r)
{
    if (r >= l) {
        int64 mid = l + (r - l) / 2;

        if (tuple3Equal(arr->data()[mid], tuple) || (tuple3Smaller(tuple, arr->data()[mid]) && tuple3Greater(tuple, arr->data()[mid-1])))
            return mid;

        if (tuple3Greater(arr->data()[mid], tuple))
            return binarySearchTuple3(arr, tuple, l, mid - 1);
 
        return binarySearchTuple3(arr, tuple, mid + 1, r);
    }
    return -1;
}


void findPivotPositions(vector<Tuple3>* arr, vector<Tuple3>* pivotsTuples, vector<int64>* pivotsPositions, int rank) {
    #pragma omp parallel for
    for (int i = 0; i < pivotsTuples->size(); i++) {
        pivotsPositions->data()[i] = binarySearchTuple3(arr, pivotsTuples->data()[i], 0, arr->size());
    }
	pivotsPositions->push_back(arr->size());
}



int getNextSendSize(int64 currentPartialPosition, int64 endPosition, int worldSize) {
    int partialSendSize = wyslijRaz; //2147483647 / worldSize;
    int64 partialSendSizeInt64 = partialSendSize;
    int64 diff = (endPosition - currentPartialPosition); 
    if (diff < partialSendSizeInt64) {
        return (int) diff;
    }
    return partialSendSize;
}


void getNextPartialPivots(vector<Tuple3>* arr, 
                          vector<Tuple3>* partialArr, 
                          vector<int64>* pivotsPosition, 
                          vector<int64>* partialPivotsPosition,
						  vector<int>* scattervPositions,
						  vector<int>* displacement,
                          int worldSize) {
    int partialArraSize = 0;
    int nextSendSize;
    
    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
        partialArraSize += nextSendSize;
    }

    partialArr->reserve(partialArraSize);
    int displacementSum = 0;

    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
		scattervPositions->push_back(nextSendSize);
        partialArr->insert(partialArr->end(), arr->begin() + partialPivotsPosition->data()[i], arr->begin() + partialPivotsPosition->data()[i] + nextSendSize);
        partialPivotsPosition->data()[i] += nextSendSize;
        displacement->data()[i] = displacementSum;
        displacementSum += nextSendSize;
    }
}


bool doNextPartialRound(vector<int64>* pivotsPosition, 
                        vector<int64>* partialPivotsPosition) {
	for (int i = 0; i < pivotsPosition->size(); i++) {
		if (pivotsPosition->data()[i] != partialPivotsPosition->data()[i]) {
			return true;
		}
	}
	return false;
}



void sendDataToProperPartition(vector<Tuple3>* A, vector<Tuple3>* A_sampleSorted, vector<int64>* pivotsPositions, int rank, int worldSize) {
    A_sampleSorted->clear();

    int nextPartitionPos = 0;
    int nextRecvNumber;

	vector<Tuple3> partialArr;
	vector<int64> partialPivotsPositions; partialPivotsPositions.resize(pivotsPositions->size());
	vector<int> scattervPositions; scattervPositions.resize(worldSize);
	vector<int> displacement; displacement.resize(worldSize);

    vector<int> arrivingNumber; arrivingNumber.resize(worldSize);
    vector<int> arrivingDisplacement; arrivingDisplacement.resize(worldSize);
    vector<Tuple3> tmp_buff;
    int sizeTmpBuff;

	partialPivotsPositions[0] = 0;
	for (int i = 1; i < pivotsPositions->size(); i++) {
		partialPivotsPositions[i] = pivotsPositions->data()[i-1];
	}

	int numberOfPartSendThisProces = ceil(((double)pivotsPositions->data()[0] / (double)wyslijRaz));
	for (int i = 1; i < pivotsPositions->size(); i++) {
		numberOfPartSendThisProces = max(numberOfPartSendThisProces, (int) ceil((double)((pivotsPositions->data()[i] - pivotsPositions->data()[i-1]) / (double)wyslijRaz)));
	}
	int currentProcesPartSends = 0;
	
    int numberOfLoops;
	// for (int p = 0; p < worldSize; p++) {

    MPI_Allreduce(&numberOfPartSendThisProces, &numberOfLoops, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);


    for (int partialSends = 0; partialSends < numberOfLoops; partialSends++) {

        partialArr.clear();
        scattervPositions.clear();
        getNextPartialPivots(A, 
                            &partialArr,
                            pivotsPositions, 
                            &partialPivotsPositions,
                            &scattervPositions,
                            &displacement,
                            worldSize);

        MPI_Alltoall((void*)scattervPositions.data(), 1, MPI_INT, (void*)arrivingNumber.data(), 1, MPI_INT, MPI_COMM_WORLD);


        sizeTmpBuff = accumulate(arrivingNumber.begin(), arrivingNumber.end(), 0);

        arrivingDisplacement.data()[0] = 0;
        for (int i = 1; i < arrivingDisplacement.size(); i++) {
            arrivingDisplacement.data()[i] = arrivingDisplacement.data()[i-1] + arrivingNumber.data()[i-1];
        }

        tmp_buff.resize(sizeTmpBuff);

        MPI_Alltoallv(partialArr.data(), 
                      scattervPositions.data(),
                      displacement.data(),
                      MPI_Tuple3,
                      tmp_buff.data(),
                      arrivingNumber.data(),
                      arrivingDisplacement.data(),
                      MPI_Tuple3,
                      MPI_COMM_WORLD);

        A_sampleSorted->insert(A_sampleSorted->end(), tmp_buff.begin(), tmp_buff.end());
        tmp_buff.clear();
    }


	
}


void sample_sort_MPI_tuple3(vector<Tuple3>* A, 
                            vector<Tuple3>* A_sampleSorted,
                            vector<Tuple3>* sample,
                            vector<Tuple3>* rootSampleRecv,
                            vector<Tuple3>* broadcastSample,
                            vector<int64>* pivotsPositions,
                            int rank, 
                            int worldSize) {

	A_sampleSorted->clear();
    sample->clear();
    broadcastSample->clear();
    broadcastSample->resize(worldSize-1);

    local_sort_openMP_tuple3(A);

    int p2 = worldSize * worldSize;

    int64 step = ceil((double) A->size() / (double) worldSize);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        sample->push_back(A->data()[minInt64(i * step, A->size()-1)]);
    }
    
    if (rank == root) {
        rootSampleRecv->resize(p2);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather((void*)sample->data(), sendNumber, MPI_Tuple3, (void*)rootSampleRecv->data(), sendNumber, MPI_Tuple3, root, MPI_COMM_WORLD);

    // // cout << "gather" << endl;s

    if (rank == root) {
        local_sort_openMP_tuple3(rootSampleRecv);
        
        for (int i = 0; i < worldSize-1; i++) {
            broadcastSample->data()[i] = rootSampleRecv->data()[(i+1) * worldSize];
        }
    }

    

    MPI_Bcast((void*)broadcastSample->data(), worldSize-1, MPI_Tuple3, root, MPI_COMM_WORLD);


    findPivotPositions(A, broadcastSample, pivotsPositions, rank);
    
	sendDataToProperPartition(A, A_sampleSorted, pivotsPositions, rank, worldSize);

	local_sort_openMP_tuple3(A_sampleSorted);
}




    // MPI_Scatter(pivotsPositions->data(), 1, MPI_Tuple3, &nextRecvNumber,
    //           num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //            MPI_Comm comm)


    // MPI_Scatterv(A->data(), 
    //              pivotsPositions->data(), 
    //              pivotsPositions->data(),
    //              MPI_Tuple3, 
    //              A_sampleSorted->data() + nextPartitionPos, 
    //              pivotsPositions->data()[rank],
    //              MPI_Tuple3, 
    //              0, 
    //              MPI_COMM_WORLD);

    

    // if (rank == 1) {
    //     cout<<
    //     for (int i = 0; i < )
    // }

    // for (int i = 0; i < worldSize; i++) {
        // MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums, num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // }

    // wyslij










// void sample_sort_MPI_Tuple3() {

// }