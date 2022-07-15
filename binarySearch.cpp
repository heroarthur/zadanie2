#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <cassert>
#include <math.h>  

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
                        int* machineWherePrefixStart,
                        int64 dataSize,
                        int64 nodeSize,
                        vector<int64>* machineSizes,
                        vector<char>* nodeCharArray,
                        int rank,
                        int worldSize) {

    int64 querySize = prefixLen;
    int64 machineWherePrefixStartLocal;

    // cout<<"preifx start get prefix "<<prefixStart<<endl;
    // cout<<"prefix start "<<prefixStart<<endl;
    int64 sizeSoFar = 0;
    for (int i = 0; i < worldSize; i++) {
        if (machineSizes->data()[i] + sizeSoFar > prefixStart) {
            machineWherePrefixStartLocal = i;
            break;
        }
        // cout<<"machine sizes "<<machineSizes->data()[i]<<endl;
        sizeSoFar += machineSizes->data()[i];
    }

    // cout<<"machine where prefix start "<<machineWherePrefixStartLocal<<endl;
    *machineWherePrefixStart = machineWherePrefixStartLocal;

    int64 offset = 0;
    for (int i = 0; i < rank; i++) {
        offset += machineSizes->data()[i];
    }


    int64 lastIndex = prefixStart + querySize + 1;
    int64 sendCount = rank < machineWherePrefixStartLocal ? 0 : minInt64(nodeSize, maxInt64(0, lastIndex - offset));
    if (rank == machineWherePrefixStartLocal) {
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
    int64 startCharArrayOffset = rank == machineWherePrefixStartLocal ? prefixStart - offset : 0;
    int sendInThisRound;
    int totalSendInThisRound;


    vector<vector<char>> prefixParts;
    if (rank == machineWherePrefixStartLocal) {
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
                    machineWherePrefixStartLocal,
                    MPI_COMM_WORLD);



        if (rank == machineWherePrefixStartLocal) {
            for (int m = 0; m < worldSize; m++) {
                prefixParts.data()[m].insert(prefixParts.data()[m].end(), receivedPrefixParts.begin() + displs.data()[m], receivedPrefixParts.begin() + displs.data()[m] + recvCountFromMachinesThisRound.data()[m]);
            }
        }

        sendCountSoFar += sendInThisRound;
    }

    if (machineWherePrefixStartLocal == rank) {
        prefix->clear();
        for (int p = 0; p < worldSize; p++) {
            prefix->insert(prefix->end(), prefixParts.data()[p].begin(), prefixParts.data()[p].end());
        }

    }


}


void getIndex(vector<int64>* machineSizes,
              vector<int64>* machineOffsets,
              int64 index,
              int* containedMachine,
              int rank,
              int worldSize) {
        
    for (int i = 0; i < worldSize; i++) {
        if (machineOffsets->data()[i] <= index && index < machineOffsets->data()[i] + machineSizes->data()[i]) {
            *containedMachine = i;
        }
    }    
}



void findAnyInfixIndexWithPrefixQuery(vector<char>* query, //ma juz \0 na koncu
                                      int64* startIndexWithPrefix,
                                      int64 dataSize,
                                      int64 nodeSize,
                                      vector<int64>* machineSizes,
                                      vector<int64>* SA_machineSizes,
                                      vector<int64>* SA_machineOffsets,
                                      vector<char>* nodeCharArray,
                                      vector<int64>* SA,
                                      int rank,
                                      int worldSize) {
    int64 l, r;
    l = 0;
    r = dataSize-1;
    int64 currIndex = (l + r) / 2;
    // cout<<"curr index "<<currIndex<<endl;
    int64 SA_index;
    int machineWhereSAstart;
    int64 prefixLen = query->size()-1;
    vector<char> prefix;
    cout<<"SA size "<<SA->size()<<endl;

    int64 sizeSoFar = 0;

    // cout<<"offsety"<<endl;
    for (int i = 0; i < worldSize; i++) {
        // cout<<SA_machineOffsets->data()[i]<<endl;
    }




    // for (int i = 0; i < worldSize; i++) {
    //     if (SA_machineOffsets->data()[i] + SA_machineSizes->data()[i] > currIndex) {
    //         machineWhereSAstart = i;
    //         break;
    //     }
    //     // cout<<"machine sizes "<<machineSizes->data()[i]<<" offset "<<SA_machineOffsets->data()[i]<<endl;
    // }



    // machineWhereSAstart = currIndex / SA_normalNodeSize;
    getIndex(SA_machineSizes,
             SA_machineOffsets,
             currIndex,
             &machineWhereSAstart,
             rank,
             worldSize);

    cout<<"machine where SA start "<<machineWhereSAstart<<" offset "<<SA_machineOffsets->data()[machineWhereSAstart]<<endl;


    // cout<<"curr index "<<currIndex<<endl;
    if (rank == machineWhereSAstart) {
        SA_index = SA->data()[currIndex - SA_machineOffsets->data()[machineWhereSAstart]];
        cout<<"wartosc sa index size "<<SA->size()<<" "<<SA_index<<endl;
    }


    MPI_Bcast(&SA_index, 1, MPI_LONG_LONG_INT, machineWhereSAstart, MPI_COMM_WORLD);
    cout<<"sa index "<<SA_index<<endl;
    // cout<<"SA index "<<SA_index<<endl;

    print_MPI_vector(SA, rank, worldSize);

    int machineWherePrefixStart;
    int cmp_dluzsze;
    int cmp_krotsze;

    for (int i = 0; i < dataSize; i++) {

        if (rank == 0) {
            // cout<<"wartosci l i r "<<l<<" "<<r<<endl;
        }
        if (l > r) {
            *startIndexWithPrefix = -1;
            return;
        }

        // if (rank == 0) {
        //     cout<<"curr index "<<currIndex<<endl;
        // }


        // cout<<"SA index "<<SA_index<<endl;
        getPrefixFromGenom(SA_index,
                           prefixLen,
                           &prefix,
                           &machineWherePrefixStart,
                           dataSize,
                           nodeSize,
                           machineSizes,
                           nodeCharArray,
                           rank,
                           worldSize);

        if (rank == 0) {
            cout<<"machine where prefix start "<<machineWherePrefixStart<<endl;
        }
        // cout<<"SAindex"<<SA_index<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (machineWherePrefixStart == rank) {
            prefix.push_back('\0');
            cmp_dluzsze = strcmp(query->data(), prefix.data());

            // cout<<"prefix dluzszy "<<prefix.data()<<" query "<<query->data()<<endl;


            prefix.pop_back(); prefix.pop_back(); prefix.push_back('\0');
            cmp_krotsze = strcmp(query->data(), prefix.data());

            // cout<<"prefix krotszy "<<prefix.data()<<"    SAindex currIndex "<<SA_index<<" "<<currIndex<<endl;
            // cout<<"cmp dluzsze "<<cmp_dluzsze<<endl;
        }

        // cout<<"machine where prefix start "<<machineWherePrefixStart<<" cmp dluzsze "<<cmp_dluzsze<<endl;

        MPI_Bcast(&cmp_dluzsze, 1, MPI_INT, machineWherePrefixStart, MPI_COMM_WORLD);
        MPI_Bcast(&cmp_krotsze, 1, MPI_INT, machineWherePrefixStart, MPI_COMM_WORLD);


        
        if (cmp_krotsze == 0 || cmp_dluzsze == 0) {
            *startIndexWithPrefix = SA_index;
            cout<<"wynik "<<currIndex<<endl;

            return;
        }



        if (cmp_dluzsze < 0) {
            r = currIndex-1;
        }
        else if (cmp_dluzsze > 0) {
            l = currIndex+1;
        }

        currIndex = (l + r) / 2;

        // if (rank == 0) {
            // cout<<"SA index przed "<<SA_index<<endl;
        // }
        // SA_index = SA->data()[currIndex];
        // if (rank == 0) {
            // cout<<"SA index po"<<SA_index<<" curr index "<<currIndex<<endl;
        // }
        // cout<<"dziwne wartosci "<<currIndex<<" "<<SA_normalNodeSize<<endl;
        // machineWhereSAstart = currIndex / SA_normalNodeSize;
        getIndex(SA_machineSizes,
                 SA_machineOffsets,
                 currIndex,
                 &machineWhereSAstart,
                 rank,
                 worldSize);

        cout<<"machine where SA start "<<machineWhereSAstart<<" offset "<<SA_machineOffsets->data()[machineWhereSAstart]<<endl;


        if (rank == machineWhereSAstart) {
            SA_index = SA->data()[currIndex - SA_machineOffsets->data()[machineWhereSAstart]];
        }
        MPI_Bcast(&SA_index, 1, MPI_LONG_LONG_INT, machineWhereSAstart, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        // break;
    }

    
}





void findMostLeftPrefix() {

}


void findMostRightPrefix() {

}



void startEdgePrefixIndexes(vector<char>* query,
                            vector<char>* nodeCharArray,
                            vector<int64>* SA,
                            int64 originalNodeSize,
                            int rank,
                            int worldSize) {

    int64 nodeSize = originalNodeSize;
    int64 nodeSAsize = SA->size();
    int64 SA_machineOffset = 0;
    int64 dataSize;

    MPI_Allreduce(&nodeSize, &dataSize, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    vector<int64> machineSizes;
    vector<int64> SA_machineSizes;
    vector<int64> SA_machineOffsets;

    machineSizes.resize(worldSize);
    MPI_Allgather(&nodeSize, 
                  1, 
                  MPI_LONG_LONG_INT, 
                  machineSizes.data(),
                  1, 
                  MPI_LONG_LONG_INT, 
                  MPI_COMM_WORLD);

    SA_machineSizes.resize(worldSize);
    MPI_Allgather(&nodeSAsize, 
                  1, 
                  MPI_LONG_LONG_INT, 
                  SA_machineSizes.data(),
                  1, 
                  MPI_LONG_LONG_INT, 
                  MPI_COMM_WORLD);

    for (int i = 0; i < rank; i++) {
        SA_machineOffset += SA_machineSizes.data()[i];
    }

    // cout<<"machine offset "<<SA_machineOffset<<endl;

    SA_machineOffsets.resize(worldSize);
    MPI_Allgather(&SA_machineOffset, 
                  1, 
                  MPI_LONG_LONG_INT, 
                  SA_machineOffsets.data(),
                  1, 
                  MPI_LONG_LONG_INT, 
                  MPI_COMM_WORLD);

    int64 startIndexWithPrefix;

    findAnyInfixIndexWithPrefixQuery(query,
                                     &startIndexWithPrefix,
                                     dataSize,
                                     nodeSize,
                                     &machineSizes,
                                     &SA_machineSizes,
                                     &SA_machineOffsets,
                                     nodeCharArray,
                                     SA,
                                     rank,
                                     worldSize);

    cout<<"znaleziono pierwsze wystapienie "<<startIndexWithPrefix<<endl;

}
    
    

                                    // (vector<char>* query, //ma juz \0 na koncu
                                    // int64* startIndexWithPrefix,
                                    // int64 dataSize,
                                    // int64 nodeSize,
                                    // vector<int64>* machineSizes,
                                    // vector<int64>* SA_machineSizes,
                                    // vector<int64>* SA_machineOffsets,
                                    // vector<char>* nodeCharArray,
                                    // vector<int64>* SA,
                                    // int rank,
                                    // int worldSize)


// CAACCTTGCGACAGGGCGGG