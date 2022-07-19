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

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include <bits/stdc++.h>

#define root 0

using namespace std;

void prepareDataForReorderSent(vector<int64>* B, 
                               vector<int64>* SA, 
                               int64 newNodeSize, 
                               int64 nodeSize,
                               int64 dataSize,
                               int64 h,
                               vector<vector<TwoInts64>>* dataForPartitions,
                               int rank,
                               int worldSize) {



        int thread_num;
        int64 threadSize;
        int nodeToSend;
        int64 offset;
        TwoInts64 data;
        vector<vector<TwoInts64>> localDataForPartitions;

        int index;
        int updateVectorsNumber;
        int vector_offset;

        //#pragma omp parallel private(thread_num, threadSize, nodeToSend, offset, data, localDataForPartitions, index, updateVectorsNumber, vector_offset)
        {           
            int threadNumber = omp_get_num_threads();

            int64 normalThreadSize = ceil(nodeSize / (double) threadNumber);
            int normalVectorUpdateNumber = ceil(worldSize / (double) threadNumber);

            thread_num = omp_get_thread_num();
            threadSize = minInt64(normalThreadSize, maxInt64(0, nodeSize - thread_num * normalThreadSize));

            localDataForPartitions.resize(worldSize);
            offset = thread_num * normalThreadSize;

            for (int64 i = offset; i < offset + threadSize; i++) {
                nodeToSend = getNodeToSend(SA->data()[i], newNodeSize);
                data.i1 = SA->data()[i];
                data.i2 = B->data()[i];
                localDataForPartitions.data()[nodeToSend].push_back(data);
            }

            index = thread_num;

            for (int i = 0; i < threadNumber; i++) {
                index = (index + i) % threadNumber;

                vector_offset = index * normalVectorUpdateNumber;
                updateVectorsNumber = minInt64(normalVectorUpdateNumber, maxInt64(0, worldSize - index * normalVectorUpdateNumber));
                for (int v_index = vector_offset; v_index < vector_offset+updateVectorsNumber; v_index++)
                {
                    dataForPartitions->data()[v_index].insert(dataForPartitions->data()[v_index].end(), 
                                                            localDataForPartitions.data()[v_index].begin(), 
                                                            localDataForPartitions.data()[v_index].end());
                }

                //#pragma omp barrier
            }
        }
}

void reorder_and_rebalance(vector<int64>* B, 
                           vector<int64>* B_new, 
                           vector<int64>* SA,
                           vector<int64>* SA_second_pointer,
                           HelpingVectorsSendingOperations* helpVectors,
                           int rank, 
                           int worldSize) {

    do_sending_operation(B, 
                         B_new, 
                         SA,
                         SA_second_pointer,
                         helpVectors,
                         EMPTY_HELP_PARAM,
                         true, 
                         rank, 
                         worldSize,
                         &prepareDataForReorderSent);
}



