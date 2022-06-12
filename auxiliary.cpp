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
#define vectorMemoryAllocationFactor 10000
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


typedef struct helpingVectorsSendingOperations {
    vector<TwoInts64> partialArr; 
    vector<int64> partialPivotsPosition;
    vector<int> scattervPositions;
    vector<int> displacement;
    vector<int> arrivingNumber;
    vector<int> arrivingDisplacement;
    vector<TwoInts64> tmp_buff;
    vector<vector<TwoInts64>> dataForPartitions; 
} HelpingVectorsSendingOperations;


typedef struct helpingVectorsSampleSort2 {
    vector<int> scattervPositions;
    vector<int> displacement;
    vector<int> arrivingNumber;
    vector<int> arrivingDisplacement;
    vector<int64> partialPivotsPosition;
    vector<int64> allArrivingNumbers;
    vector<int64> allArrivingDisplacement;
    vector<int64> addPadding;
    vector<int64> pivotsPositions;
    vector<Tuple2> partialArr; 
    vector<Tuple2> tmp_buff;
    vector<Tuple2> sample;
    vector<Tuple2> rootSampleRecv;
    vector<Tuple2> broadcastSample;
} HelpingVectorsSampleSort2;


// typedef struct helpingVectorsSampleSort3 {
//     vector<Tuple3> partialArr; 
//     vector<int64> partialPivotsPosition;
//     vector<int> scattervPositions;
//     vector<int> displacement;
//     vector<int> arrivingNumber;
//     vector<int> arrivingDisplacement;
//     vector<TwoInts64> tmp_buff; 
// } HelpingVectorsSampleSort3;



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
                          HelpingVectorsSendingOperations* helpVectors,
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
        helpVectors->dataForPartitions.data()[i].clear();
    }

    prepareDataToSent(B, SA, newNodeSize, nodeSize, dataSize, help_param, &(helpVectors->dataForPartitions), rank, worldSize);

    fill(helpVectors->partialPivotsPosition.begin(), helpVectors->partialPivotsPosition.end(), 0);

    int64 localMaxPartialSend = 0;
    int64 tmpPartialSend = 0;
    for (int i = 0; i < worldSize; i++) {
        tmpPartialSend = ceil((double) helpVectors->dataForPartitions.data()[i].size() / (double) wyslijRaz);
        localMaxPartialSend = maxInt64(localMaxPartialSend, tmpPartialSend);
    }

    int64 globalMaxPartialSend;
    
    MPI_Allreduce(&localMaxPartialSend, &globalMaxPartialSend, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);


    int sizeTmpBuff;

    for (int partialSends = 0; partialSends < globalMaxPartialSend; partialSends++) {
        getNextPartialSend(&(helpVectors->dataForPartitions), 
                           &(helpVectors->partialArr), 
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

        helpVectors->tmp_buff.resize(sizeTmpBuff);

        MPI_Alltoallv(helpVectors->partialArr.data(), 
                helpVectors->scattervPositions.data(),
                helpVectors->displacement.data(),
                MPI_TwoInts64,
                helpVectors->tmp_buff.data(),
                helpVectors->arrivingNumber.data(),
                helpVectors->arrivingDisplacement.data(),
                MPI_TwoInts64,
                MPI_COMM_WORLD);

        int64 offset = rank * newNodeSize;

        // #pragma omp parallel for
        for (int i = 0; i < helpVectors->tmp_buff.size(); i++) {
            B_help->data()[helpVectors->tmp_buff.data()[i].i1 - offset] = helpVectors->tmp_buff.data()[i].i2;
        }

        helpVectors->tmp_buff.clear();
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


void initializeHelpingVectorsSendingOperations(HelpingVectorsSendingOperations* vectors, int64 nodeDataSize, int worldSize) {
    vectors->partialArr.reserve(worldSize * wyslijRaz);
    vectors->partialPivotsPosition.resize(worldSize);
    vectors->scattervPositions.reserve(worldSize);
    vectors->displacement.reserve(worldSize);
    vectors->arrivingNumber.resize(worldSize);
    vectors->arrivingDisplacement.resize(worldSize);
    vectors->tmp_buff.reserve(worldSize * wyslijRaz);
    vectors->dataForPartitions.resize(worldSize);

    for (int i = 0; i < worldSize; i++) {
		vectors->dataForPartitions.data()[i].reserve(nodeDataSize / max(1, (worldSize-2)));
	}
}



void initializeHelpingVectorsSampleSort2(HelpingVectorsSampleSort2* vectors, int worldSize) {
    vectors->partialArr.reserve(worldSize * wyslijRaz);
    vectors->partialPivotsPosition.resize(worldSize);
    vectors->scattervPositions.reserve(worldSize);
    vectors->displacement.reserve(worldSize);
    vectors->arrivingNumber.resize(worldSize);
    vectors->arrivingDisplacement.resize(worldSize);
    vectors->allArrivingNumbers.reserve(worldSize * 2 * vectorMemoryAllocationFactor);
    vectors->allArrivingDisplacement.reserve(worldSize * 2 * vectorMemoryAllocationFactor);
    vectors->addPadding.reserve(worldSize * vectorMemoryAllocationFactor);
    vectors->tmp_buff.reserve(worldSize * wyslijRaz);
    vectors->sample.resize(worldSize);
    vectors->rootSampleRecv.resize(worldSize * worldSize);
    vectors->broadcastSample.resize(worldSize-1);
    vectors->pivotsPositions.resize(worldSize-1);   
}

