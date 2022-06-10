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

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include<bits/stdc++.h>

#define root 0


using namespace std;


void reorder_and_rebalance(vector<int64>* B, 
                           vector<int64>* B_new, 
                           vector<int64>* SA, 
                           int rank, 
                           int worldSize) {
    int64 dataSize;
    int64 nodeSize = B->size();

    MPI_Allreduce(&nodeSize, &dataSize, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    int64 newNodeSize = dataSize / worldSize;
    int64 lastNodeSize = dataSize - (worldSize-1) * newNodeSize;
    
    if (rank < worldSize-1) {
        B_new->resize(newNodeSize);
    }
    else {
        B_new->resize(lastNodeSize);
    }

    vector<vector<TwoInts64>> dataForPartitions; dataForPartitions.resize(worldSize);

    for (int i = 0; i < worldSize; i++) {
        dataForPartitions[i].clear();
    }

    #pragma omp parallel for
    for (int thread = 0; thread < worldSize; thread++) {
        for (int i = 0; i < nodeSize; i++) {
            if (thread == getNodeToSend(SA->data()[i], newNodeSize)) {
                TwoInts64 data;
                data.i1 = SA->data()[i];
                data.i2 = B->data()[i];
                dataForPartitions[thread].push_back(data);
            }
        }
        // cout<<"data size "<<dataForPartitions.data()[thread].size()<<endl;
    }

    vector<TwoInts64> partialArr; partialArr.reserve(worldSize * wyslijRaz);
    vector<int64> partialPivotsPosition; partialPivotsPosition.resize(worldSize); 
    fill(partialPivotsPosition.begin(), partialPivotsPosition.end(), 0);



    int64 localMaxPartialSend = 0;
    int64 tmpPartialSend = 0;
    for (int i = 0; i < worldSize; i++) {
        tmpPartialSend = ceil((double) dataForPartitions[i].size() / (double) wyslijRaz);
        localMaxPartialSend = maxInt64(localMaxPartialSend, tmpPartialSend);
    }


    int64 globalMaxPartialSend;
    
    MPI_Allreduce(&localMaxPartialSend, &globalMaxPartialSend, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
    // cout<<"ile sendow max "<<globalMaxPartialSend<<endl;

    vector<int> scattervPositions;
    vector<int> displacement;
    vector<int> arrivingNumber; arrivingNumber.resize(worldSize);
    vector<int> arrivingDisplacement; arrivingDisplacement.resize(worldSize);
    int sizeTmpBuff;
    vector<TwoInts64> tmp_buff; 

    for (int partialSends = 0; partialSends < globalMaxPartialSend; partialSends++) {
        getNextPartialSend(&dataForPartitions, 
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
            B_new->data()[tmp_buff.data()[i].i1 - offset] = tmp_buff.data()[i].i2;
        }

        tmp_buff.clear();
    }
}