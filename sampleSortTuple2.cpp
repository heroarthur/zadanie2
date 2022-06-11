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

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include<bits/stdc++.h>

#define root 0


using namespace std;



int64 binarySearchTuple2(vector<Tuple2>* arr, Tuple2 tuple, int64 l, int64 r)
{
    if (tuple2Greater(tuple, arr->data()[arr->size()-1])) {
        return arr->size();
    }

    if (tuple2Smaller(tuple, arr->data()[0])) {
        return 0;
    }
    // printf("tuple binary search %s   %s\n", arr->data()[arr->size()-1].B, tuple.B);


    if (r >= l) {
        int64 mid = l + (r - l) / 2;

        if (tuple2Equal(arr->data()[mid], tuple) || (tuple2Smaller(tuple, arr->data()[mid]) && tuple2Greater(tuple, arr->data()[mid-1])))
            return mid;

        if (tuple2Greater(arr->data()[mid], tuple))
            return binarySearchTuple2(arr, tuple, l, mid - 1);
 
        return binarySearchTuple2(arr, tuple, mid + 1, r);
    }
    return -1;
}


void findPivotPositionsTuple2(vector<Tuple2>* arr, vector<Tuple2>* pivotsTuples, vector<int64>* pivotsPositions, int rank) {
    // if (rank == 0)
        // cout<<"drukujemy wartosci"<<endl;

    // cout<<"arr size "<<arr->size()<<endl;
    pivotsPositions->resize(pivotsTuples->size());

    #pragma omp parallel for
    for (int i = 0; i < pivotsTuples->size(); i++) {
        if (rank == 0) {
            // printf("%s\n", pivotsTuples->data()[i].B);
        }

        pivotsPositions->data()[i] = binarySearchTuple2(arr, pivotsTuples->data()[i], 0, arr->size());
    }
	// pivotsPositions->push_back(tuple2Greater(pivotsTuples->data()[pivotsTuples->size()-1], arr->data()[arr->size()-1]) ? 0 : arr->size());
	pivotsPositions->push_back(arr->size());


    // if (rank == 3) {
    //     cout<<endl<<"PIVOTS POSITION"<<endl;
    //     for (int i = 0; i < pivotsPositions->size(); i++) {
    //         cout<<pivotsPositions->data()[i]<<" ";
    //     }
    //     cout<<endl<<endl;
    // }
}




void getNextPartialPivotsTuple2(vector<Tuple2>* arr, 
                          vector<Tuple2>* partialArr, 
                          vector<int64>* pivotsPosition, 
                          vector<int64>* partialPivotsPosition,
						  vector<int>* scattervPositions,
						  vector<int>* displacement,
                          int worldSize) {
    int partialArraSize = 0;
    int nextSendSize;

    partialArr->clear();
    scattervPositions->clear();
    
    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
        partialArraSize += nextSendSize;
    }

    partialArr->reserve(partialArraSize);
    // cout<<"partialArraySendSize "<<partialArraSize<<endl;
    int displacementSum = 0;

    for (int i = 0; i < pivotsPosition->size(); i++) {
        // cout<<"next send size "<<nextSendSize<<endl;
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
		scattervPositions->push_back(nextSendSize);
        partialArr->insert(partialArr->end(), arr->begin() + partialPivotsPosition->data()[i], arr->begin() + partialPivotsPosition->data()[i] + nextSendSize);
        partialPivotsPosition->data()[i] += nextSendSize;
        displacement->data()[i] = displacementSum;
        displacementSum += nextSendSize;
    }
}




void sendDataToProperPartitionTuple2(vector<Tuple2>* A, vector<Tuple2>* A_sampleSorted, vector<int64>* pivotsPositions, vector<int64>* allArrivingDisplacement, int rank, int worldSize) {
    A_sampleSorted->clear();
    vector<int64> allArrivingNumbers;

    int nextPartitionPos = 0;
    int nextRecvNumber;

	vector<Tuple2> partialArr;
	vector<int64> partialPivotsPositions; partialPivotsPositions.resize(pivotsPositions->size());
	vector<int> scattervPositions; scattervPositions.resize(worldSize);
	vector<int> displacement; displacement.resize(worldSize);

    vector<int> arrivingNumber; arrivingNumber.resize(worldSize);
    vector<int> arrivingDisplacement; arrivingDisplacement.resize(worldSize);
    vector<Tuple2> tmp_buff;
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

    MPI_Allreduce(&numberOfPartSendThisProces, &numberOfLoops, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);


    for (int partialSends = 0; partialSends < numberOfLoops; partialSends++) {
        // partialArr.clear();
        // scattervPositions.clear();
        getNextPartialPivotsTuple2(A, 
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
        allArrivingNumbers.insert(allArrivingNumbers.end(), arrivingNumber.begin(), arrivingNumber.end());

        tmp_buff.resize(sizeTmpBuff);

        // if (rank == 0) {
        //     cout<<"wysylamy"<<partialArr.size()<<endl;
        //     for (int i = 0; i < partialArr.size(); i++) {
        //         printf("value %s\n", partialArr.data()[i].B);
        //     }
        // }

        MPI_Alltoallv(partialArr.data(), 
                      scattervPositions.data(),
                      displacement.data(),
                      MPI_Tuple2,
                      tmp_buff.data(),
                      arrivingNumber.data(),
                      arrivingDisplacement.data(),
                      MPI_Tuple2,
                      MPI_COMM_WORLD);

        A_sampleSorted->insert(A_sampleSorted->end(), tmp_buff.begin(), tmp_buff.end());
        tmp_buff.clear();
    }

    allArrivingDisplacement->resize(allArrivingNumbers.size() + 1);
    allArrivingDisplacement->data()[0] = 0;
    for (int i = 1; i < allArrivingDisplacement->size(); i++) {
        allArrivingDisplacement->data()[i] = allArrivingDisplacement->data()[i-1] + allArrivingNumbers.data()[i-1];
    }
}

int64 roundToPowerOf2(int64 v) {
    int64 power = 1;
    while(power < v) {
        power*=2;
    }
    return power;
}


void mergeSortedParts(vector<Tuple2>* A, vector<int64>* allArrivingDisplacement, int rank) {
    int blocksNumber = allArrivingDisplacement->size()-1;
    int64 roundBlocksNumber = roundToPowerOf2(blocksNumber);
    vector<int64> addPadding; addPadding.resize(roundBlocksNumber - blocksNumber);
    fill(addPadding.begin(), addPadding.end(), allArrivingDisplacement->data()[blocksNumber]);
    allArrivingDisplacement->insert(allArrivingDisplacement->end(), addPadding.begin(), addPadding.end());
    int blocksNumberWithPadding = allArrivingDisplacement->size()-1;
	
    for (int mergeStep = 1; mergeStep < blocksNumberWithPadding; mergeStep *= 2)
	{
		int mergesInStep = (blocksNumberWithPadding / (2 * mergeStep));
        // cout<<"number of merges "<<mergesInStep<<endl;
		// #pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
            int64 indexMergeStart = 2 * mergeStep * i;
            int64 indexMergeMid = indexMergeStart + mergeStep;
            int64 indexMergeEnd = indexMergeMid + mergeStep;
            // if (rank == 0)
            //     cout<<"indeksy "<<indexMergeStart<<" "<<indexMergeMid<<" "<<indexMergeEnd<<endl;
			inplace_merge(A->begin() + allArrivingDisplacement->data()[indexMergeStart], 
                          A->begin() + allArrivingDisplacement->data()[indexMergeMid], 
                          A->begin() + allArrivingDisplacement->data()[indexMergeEnd], cmp_tuple2());
		}
	}
}



void sample_sort_MPI_tuple2(vector<Tuple2>* A, 
                            vector<Tuple2>* A_help,
                            vector<Tuple2>* sample,
                            vector<Tuple2>* rootSampleRecv,
                            vector<Tuple2>* broadcastSample,
                            vector<int64>* pivotsPositions,
                            int rank, 
                            int worldSize) {

	A_help->clear();
    sample->clear();
    broadcastSample->clear();
    broadcastSample->resize(worldSize-1);

    vector<int64> allArrivingDisplacement;

    local_sort_openMP_tuple2(A);

    int p2 = worldSize * worldSize;

    int64 step = ceil((double) A->size() / (double) worldSize);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        sample->push_back(A->data()[minInt64(i * step, A->size()-1)]);   
    }
    
    rootSampleRecv->resize(p2);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather((void*)sample->data(), sendNumber, MPI_Tuple2, (void*)rootSampleRecv->data(), sendNumber, MPI_Tuple2, MPI_COMM_WORLD);

    local_sort_openMP_tuple2(rootSampleRecv);

    for (int i = 0; i < worldSize-1; i++) {
        broadcastSample->data()[i] = rootSampleRecv->data()[(i+1) * worldSize];
    }

    findPivotPositionsTuple2(A, broadcastSample, pivotsPositions, rank);
        
	sendDataToProperPartitionTuple2(A, A_help, pivotsPositions, &allArrivingDisplacement, rank, worldSize);

    MPI_Barrier(MPI_COMM_WORLD);

    mergeSortedParts(A_help, &allArrivingDisplacement, rank);
}



