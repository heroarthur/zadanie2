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

    // #pragma omp parallel for private(data)
    // for (int thread = 0; thread < worldSize; thread++) {
        // #pragma omp parallel for private(thread, data)
        // for (int i = 0; i < nodeSize; i++) {
        //     thread = getNodeToSend(SA->data()[i], newNodeSize);
        //     // if (thread == getNodeToSend(SA->data()[i], newNodeSize)) {
        //     data.i1 = SA->data()[i];
        //     data.i2 = B->data()[i];
        //     #pragma omp critical 
        //     {
        //         dataForPartitions->data()[thread].push_back(data);
        //     }
        //     // }
        // }
        int64 normalThreadSize = ceil(nodeSize / (double) THREADS_NUM);
        // int64 lastThreadSize = maxInt64(nodeSize - (THREADS_NUM-1) * normalThreadSize, 0);

        // cout<<lastThreadSize<<endl;
    
        // cout<<"normal thread size "<<normalThreadSize<<" node size "<<nodeSize<<endl;

        #pragma omp parallel num_threads(THREADS_NUM)
        {
            int thread_num;
            int64 threadSize;
            int nodeToSend;
            int64 offset;
            TwoInts64 data;
            
            thread_num = omp_get_thread_num();
            // cout<<"thread num "<<thread_num<<endl;
            // threadSize = thread_num == THREADS_NUM-1 ? lastThreadSize : minInt64(normalThreadSize, maxInt64(0, nodeSize - (thread_num-1) * normalThreadSize));
            threadSize = minInt64(normalThreadSize, maxInt64(0, nodeSize - thread_num * normalThreadSize));

            offset = thread_num * normalThreadSize;
            // cout<<threadSize<<endl;
            vector<vector<TwoInts64>> localDataForPartitions;
            localDataForPartitions.resize(worldSize);

            for (int64 i = offset; i < offset + threadSize; i++) {
                nodeToSend = getNodeToSend(SA->data()[i], newNodeSize);
                data.i1 = SA->data()[i];
                data.i2 = B->data()[i];
                localDataForPartitions.data()[nodeToSend].push_back(data);
            }
            // cout<<"lokalny "<<localDataForPartitions.data()[0].size()<<endl;
            #pragma omp critical
            {
                for (int i = 0; i < worldSize; i++) {
                    dataForPartitions->data()[i].insert(dataForPartitions->data()[i].end(), 
                                                        localDataForPartitions.data()[i].begin(), 
                                                        localDataForPartitions.data()[i].end());
                }
            }
        }


            // #pragma omp for private(threadSize, localDataForPartitions) num_threads(THREADS_NUM)
            // for (int i = 0; i < THREADS_NUM; i++) {
            //     localDataForPartitions.resize(worldSize);

            //     int thread_num = omp_get_thread_num();
            //     nodeToSend = getNodeToSend(SA->data()[i], newNodeSize);
            //     data.i1 = SA->data()[i];
            //     data.i2 = B->data()[i];
            //     localDataForPartitions.data()[nodeToSend].push_back(data);
            // }

            // #pragma omp critical
            // {
            //     for (int i = 0; i < worldSize; i++) {
            //         dataForPartitions->data()[i].insert(dataForPartitions->data()[i].end(), localDataForPartitions.data()[i].begin(), localDataForPartitions.data()[i].end());
            //     }
            // }

        // #pragma omp parallel
        // {
        //     int tmpValue = 0;
        //     #pragma omp for nowait
        //     for (int i = 0; i < inputArray.size(); i++) {
        //     tmpValue += inputArray[i];
        //     }
        //     #pragma omp atomic
        //     totalValue += tmpValue;
        // }
    // }
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



