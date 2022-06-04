#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>
#include <mpi.h>

#ifndef   sortsOpenMP
#define   sortsOpenMP

#include "sortsOpenMP.cpp"

#endif

#include<bits/stdc++.h>

#define root 0

#define wyslijRaz 5

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

bool tuple3Equal(Tuple3 t1, Tuple3 t2) {
    return t1.B == t2.B && t1.B2 == t2.B2;
}

bool tuple3Greater(Tuple3 t1, Tuple3 t2) {
    return t1.B > t2.B || (t1.B == t2.B && t1.B2 > t2.B2);
}

bool tuple3Smaller(Tuple3 t1, Tuple3 t2) {
    return !tuple3Equal(t1, t2) && !tuple3Greater(t1, t2);
}


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


void findPivotPositions(vector<Tuple3>* arr, vector<Tuple3>* pivotsTuples, vector<int64>* pivotsPositions) {
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
	displacement->data()[0] = 0;


    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
		scattervPositions->push_back(nextSendSize);
        partialArr->insert(partialArr->end(), arr->begin() + partialPivotsPosition->data()[i], arr->begin() + partialPivotsPosition->data()[i] + nextSendSize);
        partialPivotsPosition->data()[i] += nextSendSize;
    }

	for (int i = 1; i < pivotsPosition->size(); i++) {
		displacement->data()[i] = displacement->data()[i-1] + scattervPositions->data()[i];
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

    int nextPartitionPos = 0;
    int nextRecvNumber;

	vector<Tuple3> partialArr;
	vector<int64> partialPivotsPositions; partialPivotsPositions.resize(pivotsPositions->size());
	vector<int> scattervPositions; scattervPositions.resize(worldSize);
	vector<int> displacement; displacement.resize(worldSize);
	partialPivotsPositions[0] = 0;
	for (int i = 1; i < pivotsPositions->size(); i++) {
		partialPivotsPositions[i] = pivotsPositions->data()[i-1];
	}

	int numberOfPartSendThisProces = ceil(((double)pivotsPositions->data()[0] / (double)wyslijRaz));
	for (int i = 1; i < pivotsPositions->size(); i++) {
		numberOfPartSendThisProces = max(numberOfPartSendThisProces, (int) ceil((double)((pivotsPositions->data()[i] - pivotsPositions->data()[i-1]) / (double)wyslijRaz)));
	}

	int currentProcesPartSends = 0;

	
	for (int p = 0; p < 1; p++) {

		if (rank == p) {
			currentProcesPartSends = numberOfPartSendThisProces;
		}

		MPI_Bcast(&currentProcesPartSends, 1, MPI_INT, p, MPI_COMM_WORLD);

		for (int sends = 0; sends < currentProcesPartSends; sends++) {

			if (rank == 0) {
			partialArr.clear();
			scattervPositions.clear();
			getNextPartialPivots(A, 
								&partialArr,
								pivotsPositions, 
								&partialPivotsPositions,
								&scattervPositions,
								&displacement,
								worldSize);


			}

			MPI_Scatter(scattervPositions.data(), 1, MPI_INT, &nextRecvNumber, 1, MPI_INT, p, MPI_COMM_WORLD);

			int64 previousSize = A_sampleSorted->size();
			A_sampleSorted->resize(previousSize + nextRecvNumber);

			MPI_Scatterv(partialArr.data(), 
						scattervPositions.data(), 
						displacement.data(),
						MPI_Tuple3, 
						A_sampleSorted->data() + previousSize, 
						nextRecvNumber,
						MPI_Tuple3, 
						0, 
						MPI_COMM_WORLD);
		}
	}
}


void sample_sort_MPI_tuple3(vector<Tuple3>* A, 
                            vector<Tuple3>* A_sampleSorted,
                            vector<Tuple3>* sample,
                            vector<Tuple3>* rootSampleRecv,
                            vector<Tuple3>* broadcastSample,
                            vector<int64>* pivotsPositions,
                            int64* sizeOnNode, 
                            int rank, 
                            int worldSize) {

    local_sort_openMP_tuple3(A);
    // cout<<"local sorted "<<worldSize<<" "<<rank<<endl;

	A_sampleSorted->clear();

    int p2 = worldSize * worldSize;


    int64 step = *sizeOnNode / (worldSize + 1);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        sample->data()[i] = A->data()[(i + 1) * step];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(sample->data(), sendNumber, MPI_Tuple3, (void*)rootSampleRecv->data(), sendNumber, MPI_Tuple3, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // // cout << "gather" << endl;s

    if (rank == root) {
        local_sort_openMP_tuple3(rootSampleRecv);
        for (int i = 0; i < worldSize-1; i++) {
            broadcastSample->data()[i] = rootSampleRecv->data()[(i+1) * worldSize];
        }
    }

    

    MPI_Bcast((void*)broadcastSample->data(), worldSize-1, MPI_Tuple3, root, MPI_COMM_WORLD);


    findPivotPositions(A, broadcastSample, pivotsPositions);
    
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