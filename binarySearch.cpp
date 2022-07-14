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






void getPrefixFromGenom(int64 prefixStart,
                        int64 prefixLen,
                        vector<char>* prefix,
                        int64 dataSize,
                        int64 nodeSize,
                        vector<int64>* machineSizes,
                        int64 newNodeSize,
                        vector<char>* nodeCharArray,
                        int rank,
                        int worldSize) {

    int64 querySize = prefixLen;


    int machineWherePrefixStart;
    int64 sizeSoFar = 0;
    for (int i = 0; i < worldSize; i++) {
        if (machineSizes->data()[i] + sizeSoFar > prefixStart) {
            machineWherePrefixStart = i;
            break;
        }
        sizeSoFar += machineSizes->data()[i];
    }

    int64 offset = 0;
    for (int i = 0; i < rank; i++) {
        offset += machineSizes->data()[i];
    }


    int64 lastIndex = prefixStart + querySize + 1;
    int64 sendCount = rank < machineWherePrefixStart ? 0 : minInt64(nodeSize, maxInt64(0, lastIndex - offset));
    if (rank == machineWherePrefixStart) {
        sendCount = minInt64(lastIndex - prefixStart, nodeSize - (prefixStart - offset));
    }

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
    int64 startCharArrayOffset = rank == machineWherePrefixStart ? prefixStart - offset : 0;
    int sendInThisRound;
    int totalSendInThisRound;


    vector<vector<char>> prefixParts;
    if (rank == machineWherePrefixStart) {
        prefixParts.resize(worldSize);
    }

    vector<char> receivedPrefixParts;


    
    for (int i = 0; i < globalMaxPartialSend; i++) {

        sendInThisRound = minInt64(wyslijRaz, sendCount - sendCountSoFar);

        MPI_Allgather(&sendInThisRound, 
                      1, 
                      MPI_INT, 
                      recvCountFromMachinesThisRound.data(),
                      1, 
                      MPI_INT, 
                      MPI_COMM_WORLD);

        displs.data()[0] = 0;
        for (int d = 1; d < worldSize; d++) {
            displs.data()[d] = displs.data()[d-1] + recvCountFromMachinesThisRound.data()[d-1];
        }

        totalSendInThisRound = 0;
        MPI_Allreduce(&sendInThisRound, &totalSendInThisRound, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        receivedPrefixParts.resize(totalSendInThisRound);


        MPI_Gatherv(nodeCharArray->data() + sendCountSoFar + startCharArrayOffset,
                    sendInThisRound, 
                    MPI_CHAR, 
                    receivedPrefixParts.data(),
                    recvCountFromMachinesThisRound.data(), 
                    displs.data(), 
                    MPI_CHAR, 
                    machineWherePrefixStart,
                    MPI_COMM_WORLD);



        if (rank == machineWherePrefixStart) {
            for (int m = 0; m < worldSize; m++) {
                prefixParts.data()[m].insert(prefixParts.data()[m].end(), receivedPrefixParts.begin() + displs.data()[m], receivedPrefixParts.begin() + displs.data()[m] + recvCountFromMachinesThisRound.data()[m]);
            }
        }

        sendCountSoFar += sendInThisRound;
    }

    if (machineWherePrefixStart == rank) {
        prefix->clear();
        for (int p = 0; p < worldSize; p++) {
            prefix->insert(prefix->end(), prefixParts.data()[p].begin(), prefixParts.data()[p].end());
        }

    }


}




