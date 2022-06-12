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
#define charArrayLen 8
#define EMPTY_HELP_PARAM 0
 //(2147483647 / 10)

using namespace std;


typedef struct mpi_tuple2 {
    char B[charArrayLen+1];
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


typedef struct helpingVectors {
    vector<TwoInts64> partialArr; 
    vector<int64> partialPivotsPosition;
    vector<int> scattervPositions;
    vector<int> displacement;
    vector<int> arrivingNumber;
    vector<int> arrivingDisplacement;
    vector<TwoInts64> tmp_buff; 
} HelpingVectors;


struct cmp_tuple3 {
    bool operator ()(Tuple3 const& a, Tuple3 const& b) const {
        return a.B < b.B || (a.B == b.B && a.B2 < b.B2);
    }
};


struct cmp_tuple2 {
    bool operator ()(Tuple2 const& a, Tuple2 const& b) const {
        return strcmp(a.B, b.B) < 0;
    }
};

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
    // cout<<"end position "<<endPosition<<endl;
    int64 diff = (endPosition - currentPartialPosition); 
    if (diff < partialSendSizeInt64) {
        return (int) diff;
    }
    return partialSendSize;
}



void initialize_SA(vector<int64>* __restrict__ SA, 
                   vector<Tuple2>* __restrict__ tuple2) {
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

inline int getNodeToSend(int64 id, int64 nodeSize) {
    return id / nodeSize;
}




void do_sending_operation(vector<int64>* B, 
                          vector<int64>* B_help, 
                          vector<int64>* SA,
                          vector<vector<TwoInts64>>* dataForPartitions,
                          HelpingVectors* helpVectors,
                          int64 help_param, 
                          int rank, 
                          int worldSize,
                          void (*prepareDataToSent)(vector<int64>*, vector<int64>*, int64, int64, int64, int64, vector<vector<TwoInts64>>*, int, int)) {
    int64 dataSize;
    int64 nodeSize = B->size();

    MPI_Allreduce(&nodeSize, &dataSize, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    int64 newNodeSize = dataSize / worldSize;
    int64 lastNodeSize = dataSize - (worldSize-1) * newNodeSize;
    
    if (rank < worldSize-1) {
        B_help->resize(newNodeSize);
    }
    else {
        B_help->resize(lastNodeSize);
    }

    for (int i = 0; i < worldSize; i++) {
        dataForPartitions->data()[i].clear();
    }

    prepareDataToSent(B, SA, newNodeSize, nodeSize, dataSize, help_param, dataForPartitions, rank, worldSize);

    vector<TwoInts64> partialArr; partialArr.reserve(worldSize * wyslijRaz);
    vector<int64> partialPivotsPosition; partialPivotsPosition.resize(worldSize); 
    fill(partialPivotsPosition.begin(), partialPivotsPosition.end(), 0);

    int64 localMaxPartialSend = 0;
    int64 tmpPartialSend = 0;
    for (int i = 0; i < worldSize; i++) {
        tmpPartialSend = ceil((double) dataForPartitions->data()[i].size() / (double) wyslijRaz);
        localMaxPartialSend = maxInt64(localMaxPartialSend, tmpPartialSend);
    }

    int64 globalMaxPartialSend;
    
    MPI_Allreduce(&localMaxPartialSend, &globalMaxPartialSend, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);

    vector<int> scattervPositions;
    vector<int> displacement;
    vector<int> arrivingNumber; arrivingNumber.resize(worldSize);
    vector<int> arrivingDisplacement; arrivingDisplacement.resize(worldSize);
    int sizeTmpBuff;
    vector<TwoInts64> tmp_buff; 

    for (int partialSends = 0; partialSends < globalMaxPartialSend; partialSends++) {
        getNextPartialSend(dataForPartitions, 
                           &partialArr, 
                           &partialPivotsPosition,
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
                MPI_TwoInts64,
                tmp_buff.data(),
                arrivingNumber.data(),
                arrivingDisplacement.data(),
                MPI_TwoInts64,
                MPI_COMM_WORLD);

        int64 offset = rank * newNodeSize;

        // #pragma omp parallel for
        for (int i = 0; i < tmp_buff.size(); i++) {
            B_help->data()[tmp_buff.data()[i].i1 - offset] = tmp_buff.data()[i].i2;
        }

        tmp_buff.clear();
    }
}


void print_MPI_vector(vector<int64>* v, int rank, int worldSize) { 
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < worldSize; r++) {
        if (r == rank) {
		    for (int i = 0; i < v->size(); i++) {
		    	printf("%lld ", v->data()[i]);
		    }
	    }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        cout<<endl<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
}



void print_MPI_tuple2(vector<Tuple2>* v, int rank, int worldSize) { 
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < worldSize; r++) {
        if (r == rank) {
		    for (int i = 0; i < v->size(); i++) {
		    	printf("%s\n", v->data()[i].B);
		    }
	    }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        cout<<endl<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
}


void switchPointersTuple2(vector<Tuple2>** A1, vector<Tuple2>** A2) {
    vector<Tuple2> *A_tmp_pointer;
    A_tmp_pointer = *A2;
	*A2 = *A1;
	*A1 = A_tmp_pointer;
}

void switchPointersInt64(vector<int64>** A1, vector<int64>** A2) {
    vector<int64> *A_tmp_pointer;
    A_tmp_pointer = *A2;
	*A2 = *A1;
	*A1 = A_tmp_pointer;
}


void initializeHelpingVectors(HelpingVectors* vectors, int worldSize) {
    vectors->partialArr.reserve(worldSize * wyslijRaz);
    vectors->partialPivotsPosition.resize(worldSize);
    vectors->scattervPositions.reserve(worldSize);
    vectors->displacement.reserve(worldSize);
    vectors->arrivingNumber.resize(worldSize);
    vectors->arrivingDisplacement.resize(worldSize);
    vectors->tmp_buff.reserve(worldSize * wyslijRaz);
}
