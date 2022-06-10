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


#include <bits/stdc++.h>

#define root 0
#define wyslijRaz (2147483647 / 10)
#define int64 long long int
#define charArrayLen 4
 //(2147483647 / 10)

using namespace std;


typedef struct mpi_tuple2 {
    char B[charArrayLen];
    int64 i;
} Tuple2;


typedef struct mpi_tuple3 {
    int64 B;
    int64 B2;
    int64 i;
} Tuple3;


typedef struct twoInts64 {
    int64 i1;
    int64 i2;
} TwoInts64;

MPI_Datatype MPI_Tuple2;
MPI_Datatype MPI_Tuple3;

MPI_Datatype MPI_TwoInts64;


int get_block_start(int blockId, int blocksNumber, int dataSize) {
	return (dataSize / blocksNumber) * blockId;
}

bool tuple3Equal(Tuple3 t1, Tuple3 t2) {
    return t1.B == t2.B && t1.B2 == t2.B2;
}

bool tuple3Greater(Tuple3 t1, Tuple3 t2) {
    return t1.B > t2.B || (t1.B == t2.B && t1.B2 > t2.B2);
}

bool tuple3Smaller(Tuple3 t1, Tuple3 t2) {
    return !tuple3Equal(t1, t2) && !tuple3Greater(t1, t2);
}


bool tuple2Equal(Tuple2 t1, Tuple2 t2) {
	return strcmp(t1.B, t2.B) == 0;
}

bool tuple2Greater(Tuple2 t1, Tuple2 t2) {
    return strcmp(t1.B, t2.B) > 0;
}

bool tuple2Smaller(Tuple2 t1, Tuple2 t2) {
    return !tuple2Equal(t1, t2) && !tuple2Greater(t1, t2);
}


int64 minInt64(int64 a, int64 b) {
	return a < b ? a : b;
}

int64 maxInt64(int64 a, int64 b) {
	return a > b ? a : b;
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


int getNextSendSize(int64 currentPartialPosition, int64 endPosition, int worldSize) {
    int partialSendSize = wyslijRaz; //2147483647 / worldSize;
    int64 partialSendSizeInt64 = partialSendSize;
    int64 diff = (endPosition - currentPartialPosition); 
    if (diff < partialSendSizeInt64) {
        return (int) diff;
    }
    return partialSendSize;
}



void initialize_SA(vector<int64>* SA, vector<Tuple2>* tuple2) {
    SA->resize(tuple2->size());

    #pragma omp parallel for
    for (int i = 0; i < tuple2->size(); i++) {
        SA->data()[i] = tuple2->data()[i].i;
    }
}




void getNextPartialSend(vector<vector<TwoInts64>>* dataForPartitions, 
                        vector<TwoInts64>* partialArr, 
                        vector<int64>* partialPivotsPosition,
						vector<int>* scattervPositions,
						vector<int>* displacement,
                        int worldSize) {

    partialArr->clear();
    scattervPositions->clear();
    displacement->clear();

    int pivot = 0;
    
    for(int node = 0; node < worldSize; node++) {
        pivot = 0;
        // cout<<"data for partition "<<dataForPartitions->data()[node].size()<<endl;
        for(int i = partialPivotsPosition->data()[node]; i < minInt64(partialPivotsPosition->data()[node] + wyslijRaz, dataForPartitions->data()[node].size()); i++) {
            partialArr->push_back(dataForPartitions->data()[node].data()[i]);
            pivot++;
        }
        // cout<<"pivot size "<<pivot<<endl;
        scattervPositions->push_back(pivot);
    }
    displacement->push_back(0);
    int nextDispl;
    for(int i = 1; i < worldSize; i++) {
        nextDispl = displacement->data()[i-1] + scattervPositions->data()[i-1];
        displacement->push_back(nextDispl);
    }
}

bool doNextPartialSend(vector<int64>* pivotsPosition, 
                       vector<int64>* partialPivotsPosition) {
    
    for(int i = 0; i < pivotsPosition->size(); i++) {
        if (pivotsPosition->data()[i] > partialPivotsPosition->data()[i]) {
            return true;
        }
    }
    return false;
}

int getNodeToSend(int64 id, int64 nodeSize) {
    return id / nodeSize;
}