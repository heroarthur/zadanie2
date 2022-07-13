#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <cassert>

#include <omp.h>
#include <mpi.h>

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif



    
// int64 dataSize;
// int64 nodeSize = A->size();

// MPI_Allreduce(&nodeSize, &dataSize, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

// int64 newNodeSize = ceil(dataSize / (double) worldSize);
// int64 thisNodeNewSize = minInt64(newNodeSize, maxInt64(0, dataSize - rank * newNodeSize));






void getPrefixFromGenom(vector<char>* query, 
                        vector<char>* prefix,
                        int64 dataSize,
                        int64 nodeSize,
                        int64 newNodeSize,
                        vector<int64> machinesSizes,
                        int64 prefixStart,
                        vector<int64>* nodeCharArray,
                        HelpingVectorsSendingOperations* helpVectors,
                        int rank,
                        int worldSize) {

    int64 querySize = query->size();

    int machineWherePrefixStart = prefixStart / newNodeSize;
    int64 offset = rank * newNodeSize;

    if (machineWherePrefixStart == rank) {
        prefix->resize(querySize);
    }
    else {
        prefix->clear();
    }


    int64 lastIndex = prefixStart + querySize + 1;
    int64 sendCount = rank < machineWherePrefixStart ? 0 : minInt64(nodeSize, lastIndex - offset);
    vector<int64> totalRecvCountFromMachines;
    totalRecvCountFromMachines.resize(worldSize);
    vector<int> recvCountFromMachinesThisRound;
    recvCountFromMachinesThisRound.resize(worldSize);

    vector<int> displs;
    displs.resize(worldSize);

    int64 localMaxPartialSend = ceil((double) sendCount / (double) wyslijRaz);

    int64 globalMaxPartialSend;
    
    MPI_Allreduce(&localMaxPartialSend, &globalMaxPartialSend, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);

    MPI_Allgather(&sendCount, 
                  1, 
                  MPI_LONG_LONG_INT, 
                  totalRecvCountFromMachines.data(),
                  1, 
                  MPI_LONG_LONG_INT, 
                  MPI_COMM_WORLD);


    int64 sendCountSoFar = 0;
    int64 sendInThisRound;
    int64 charrArrStartOffset = rank > machineWherePrefixStart ? 0 : prefixStart - offset;
    
    for (int i = 0; i < globalMaxPartialSend; i++) {

        sendInThisRound = maxInt64(wyslijRaz, sendCount - sendCountSoFar); 
        
        MPI_Allgather(&sendInThisRound, 
                      1, 
                      MPI_LONG_LONG_INT, 
                      recvCountFromMachinesThisRound.data(),
                      1, 
                      MPI_LONG_LONG_INT, 
                      MPI_COMM_WORLD);

        displs.data()[0] = 0;
        for (int d = 1; d < worldSize; d++) {
            displs.data()[d] = displs.data()[d-1] + recvCountFromMachinesThisRound.data()[d-1] + sendCountSoFar;
        }

        MPI_Gatherv(nodeCharArray->data() + charrArrStartOffset,
                    sendCount, 
                    MPI_CHAR, 
                    prefix->data(),
                    recvCountFromMachinesThisRound.data(), 
                    displs.data(), 
                    MPI_CHAR, 
                    machineWherePrefixStart,
                    MPI_COMM_WORLD);

        sendCountSoFar += sendInThisRound;
    }
}


// int64 findAnyIndexWithEqualPrefix(vector<char>* query,
//                                   int64* equalPrefixIndex,
//                                   int64 currIndex) {
    
//     // zdobadz prefix z tej pozycji

//     // jezeli jest to zwroc jak nie to szukaj dalej


// }


