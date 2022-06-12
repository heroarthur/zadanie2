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



int64 binarySearchTuple2(vector<Tuple2>* arr, 
                         Tuple2 tuple, 
                         int64 l, 
                         int64 r)
{
    if (tuple2Greater(tuple, arr->data()[arr->size()-1])) {
        return arr->size();
    }

    if (tuple2Smaller(tuple, arr->data()[0])) {
        return 0;
    }

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


void findPivotPositionsTuple2(vector<Tuple2>* arr, 
                              vector<Tuple2>* pivotsTuples, 
                              vector<int64>* pivotsPositions, 
                              int rank) {

    pivotsPositions->resize(pivotsTuples->size());

    #pragma omp parallel for
    for (int i = 0; i < pivotsTuples->size(); i++) {
        pivotsPositions->data()[i] = binarySearchTuple2(arr, pivotsTuples->data()[i], 0, arr->size());
    }
	pivotsPositions->push_back(arr->size());
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
    
    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
        partialArraSize += nextSendSize;
    }

    int displacementSum = 0;

    for (int i = 0; i < pivotsPosition->size(); i++) {
        nextSendSize = getNextSendSize(partialPivotsPosition->data()[i], pivotsPosition->data()[i], worldSize);
		scattervPositions->data()[i] = nextSendSize;
        partialArr->insert(partialArr->end(), arr->begin() + partialPivotsPosition->data()[i], arr->begin() + partialPivotsPosition->data()[i] + nextSendSize);
        partialPivotsPosition->data()[i] += nextSendSize;
        displacement->data()[i] = displacementSum;
        displacementSum += nextSendSize;
    }
}



void sendDataToProperPartitionTuple2(vector<Tuple2>* A, 
                                     vector<Tuple2>* A_sampleSorted, 
                                     HelpingVectorsSampleSort2* helpVectors,
                                     int rank, 
                                     int worldSize) {

    A_sampleSorted->clear();

    int nextPartitionPos = 0;
    int nextRecvNumber;

	helpVectors->partialPivotsPosition.resize(helpVectors->pivotsPositions.size());
	helpVectors->scattervPositions.resize(worldSize);
	helpVectors->displacement.resize(worldSize);
    helpVectors->arrivingNumber.resize(worldSize);
    helpVectors->arrivingDisplacement.resize(worldSize);
    
    int sizeTmpBuff;

	helpVectors->partialPivotsPosition[0] = 0;
	for (int i = 1; i < helpVectors->pivotsPositions.size(); i++) {
		helpVectors->partialPivotsPosition[i] = helpVectors->pivotsPositions.data()[i-1];
	}

	int numberOfPartSendThisProces = ceil(((double)helpVectors->pivotsPositions.data()[0] / (double)wyslijRaz));
	for (int i = 1; i < helpVectors->pivotsPositions.size(); i++) {
		numberOfPartSendThisProces = max(numberOfPartSendThisProces, (int) ceil((double)((helpVectors->pivotsPositions.data()[i] - helpVectors->pivotsPositions.data()[i-1]) / (double)wyslijRaz)));
	}
	int currentProcesPartSends = 0;
	
    int numberOfLoops;

    MPI_Allreduce(&numberOfPartSendThisProces, &numberOfLoops, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);


    for (int partialSends = 0; partialSends < numberOfLoops; partialSends++) {

        getNextPartialPivotsTuple2(A, 
                            &(helpVectors->partialArr),
                            &(helpVectors->pivotsPositions), 
                            &(helpVectors->partialPivotsPosition),
                            &(helpVectors->scattervPositions),
                            &(helpVectors->displacement),
                            worldSize);

        MPI_Alltoall((void*)helpVectors->scattervPositions.data(), 1, MPI_INT, (void*)helpVectors->arrivingNumber.data(), 1, MPI_INT, MPI_COMM_WORLD);


        sizeTmpBuff = accumulate(helpVectors->arrivingNumber.begin(), helpVectors->arrivingNumber.end(), 0);

        helpVectors->arrivingDisplacement.data()[0] = 0;
        for (int i = 1; i < helpVectors->arrivingDisplacement.size(); i++) {
            helpVectors->arrivingDisplacement.data()[i] = helpVectors->arrivingDisplacement.data()[i-1] + helpVectors->arrivingNumber.data()[i-1];
        }
        helpVectors->allArrivingNumbers.insert(helpVectors->allArrivingNumbers.end(), helpVectors->arrivingNumber.begin(), helpVectors->arrivingNumber.end());

        helpVectors->tmp_buff.resize(sizeTmpBuff);

        MPI_Alltoallv(helpVectors->partialArr.data(), 
                      helpVectors->scattervPositions.data(),
                      helpVectors->displacement.data(),
                      MPI_Tuple2,
                      helpVectors->tmp_buff.data(),
                      helpVectors->arrivingNumber.data(),
                      helpVectors->arrivingDisplacement.data(),
                      MPI_Tuple2,
                      MPI_COMM_WORLD);

        A_sampleSorted->insert(A_sampleSorted->end(), helpVectors->tmp_buff.begin(), helpVectors->tmp_buff.end());
        helpVectors->tmp_buff.clear();
    }

    helpVectors->allArrivingDisplacement.resize(helpVectors->allArrivingNumbers.size() + 1);
    helpVectors->allArrivingDisplacement.data()[0] = 0;
    for (int i = 1; i < helpVectors->allArrivingDisplacement.size(); i++) {
        helpVectors->allArrivingDisplacement.data()[i] = helpVectors->allArrivingDisplacement.data()[i-1] + helpVectors->allArrivingNumbers.data()[i-1];
    }
}

int64 roundToPowerOf2(int64 v) {
    int64 power = 1;
    while(power < v) {
        power*=2;
    }
    return power;
}


void mergeSortedParts(vector<Tuple2>* A, 
                      HelpingVectorsSampleSort2* helpVectors, 
                      int rank) {
                          
    int blocksNumber = helpVectors->allArrivingDisplacement.size()-1;
    int64 roundBlocksNumber = roundToPowerOf2(blocksNumber);
    helpVectors->addPadding.resize(roundBlocksNumber - blocksNumber);
    fill(helpVectors->addPadding.begin(), helpVectors->addPadding.end(), helpVectors->allArrivingDisplacement.data()[blocksNumber]);
    helpVectors->allArrivingDisplacement.insert(helpVectors->allArrivingDisplacement.end(), helpVectors->addPadding.begin(), helpVectors->addPadding.end());
    int blocksNumberWithPadding = helpVectors->allArrivingDisplacement.size()-1;

    for (int mergeStep = 1; mergeStep < blocksNumberWithPadding; mergeStep *= 2)
	{
		int mergesInStep = (blocksNumberWithPadding / (2 * mergeStep));

		for (int i = 0; i < mergesInStep; i++) {
            int64 indexMergeStart = 2 * mergeStep * i;
            int64 indexMergeMid = indexMergeStart + mergeStep;
            int64 indexMergeEnd = indexMergeMid + mergeStep;



			inplace_merge(A->begin() + helpVectors->allArrivingDisplacement.data()[indexMergeStart], 
                          A->begin() + helpVectors->allArrivingDisplacement.data()[indexMergeMid], 
                          A->begin() + helpVectors->allArrivingDisplacement.data()[indexMergeEnd], cmp_tuple2());
		}
	}
}



void sample_sort_MPI_tuple2(vector<Tuple2>* A, 
                            vector<Tuple2>* A_help,
                            HelpingVectorsSampleSort2* helpVectors,
                            int rank, 
                            int worldSize) {

	A_help->clear();
    helpVectors->allArrivingNumbers.clear();

    local_sort_openMP_tuple2(A);

    int p2 = worldSize * worldSize;

    int64 step = ceil((double) A->size() / (double) worldSize);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        helpVectors->sample.data()[i] = A->data()[minInt64(i * step, A->size()-1)];   
    }
        
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather((void*)helpVectors->sample.data(), sendNumber, MPI_Tuple2, (void*)helpVectors->rootSampleRecv.data(), sendNumber, MPI_Tuple2, MPI_COMM_WORLD);

    local_sort_openMP_tuple2(&(helpVectors->rootSampleRecv));

    for (int i = 0; i < worldSize-1; i++) {
        helpVectors->broadcastSample.data()[i] = helpVectors->rootSampleRecv.data()[(i+1) * worldSize];
    }

    findPivotPositionsTuple2(A, &(helpVectors->broadcastSample), &(helpVectors->pivotsPositions), rank);
        
	sendDataToProperPartitionTuple2(A, A_help, helpVectors, rank, worldSize);

    MPI_Barrier(MPI_COMM_WORLD);

    mergeSortedParts(A_help, helpVectors, rank);
}



