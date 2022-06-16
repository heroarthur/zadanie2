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


#include<bits/stdc++.h>

#define root 0


using namespace std;



void rebucket_assign_2h_group_rank(vector<Tuple3>* __restrict__ tuple, 
                               vector<int64>* __restrict__ B,
                               vector<int64>* __restrict__ SA, 
                               bool* allSingletones,
                               int rank,
                               int worldSize) {

    bool allSingletonesCheck;
    int64 indexOffset;
    Tuple3 lastTuple;

    if (rank == 0) {
        bool localAllSingletones = true;
        int64 size = tuple->size();
        int64 lastPossibleIndex = size-1;
        B->resize(size);
        SA->resize(size);
        B->data()[0] = 0;
        SA->data()[0] = tuple->data()[0].i;
        for (int64 i = 1; i < size; i++) {
            if (tuple3Equal(tuple->data()[i-1], tuple->data()[i])) {
                    localAllSingletones = false;
                    B->data()[i] = B->data()[i-1];
                }
            else {
                B->data()[i] = i;
            }
            SA->data()[i] = tuple->data()[i].i;
        }

        *allSingletones = localAllSingletones;
        indexOffset = lastPossibleIndex + 1;
        lastTuple = tuple->data()[lastPossibleIndex];
        
        if (worldSize > 1) {
            MPI_Send(allSingletones, 1, MPI_C_BOOL, 1, root, MPI_COMM_WORLD);
            MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, 1, root, MPI_COMM_WORLD);
            MPI_Send(&lastTuple, 1, MPI_Tuple3, 1, root, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(allSingletones, 1, MPI_C_BOOL,        rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&indexOffset, 1,         MPI_LONG_LONG_INT, rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&lastTuple, 1,           MPI_Tuple3,        rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int64 size = tuple->size();
        int64 lastPossibleIndex = size-1;
        B->resize(size);
        SA->resize(size);
        B->data()[0] = tuple3Equal(lastTuple, tuple->data()[0]) ? indexOffset-1 : indexOffset;
        SA->data()[0] = tuple->data()[0].i;
        for (int64 i = 1; i < size; i++) {
            cout<<"B ("<<tuple->data()[i-1].B<<","<<tuple->data()[i-1].B2<<") B2 ("<<tuple->data()[i].B<<", "<<tuple->data()[i].B2<<")"<<endl;

            if (tuple3Equal(tuple->data()[i-1], tuple->data()[i])) {
                    *allSingletones = false;
                    B->data()[i] = B->data()[i-1];
                    // break;
                }
            else {
                B->data()[i] = i + indexOffset;
            }
            SA->data()[i] = tuple->data()[i].i;

        }

        indexOffset += lastPossibleIndex + 1;
        lastTuple = tuple->data()[lastPossibleIndex];

        if (rank < worldSize-1) {
            MPI_Send(allSingletones, 1, MPI_C_BOOL, rank+1, rank, MPI_COMM_WORLD);
            MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, rank+1, rank, MPI_COMM_WORLD);
            MPI_Send(&lastTuple, 1, MPI_Tuple3, rank+1, rank, MPI_COMM_WORLD);
        }
        else 
        {
            cout<<"czy wszystkie singletone "<<*allSingletones<<endl;
        }
    }

    MPI_Bcast(allSingletones, 1, MPI_C_BOOL, worldSize-1, MPI_COMM_WORLD);
}



void rebucket_assign_h_group_rank(vector<Tuple2>* __restrict__ tuple, 
                                  vector<int64>* __restrict__ B, 
                                  int rank,
                                  int worldSize) {

    int64 indexOffset;
    Tuple2 lastTuple;

    if (rank == 0) {
        int64 size = tuple->size();
        int64 lastPossibleIndex = size-1;
        B->resize(size);
        B->data()[0] = 0;

        for (int64 i = 1; i < size; i++) {
            if (tuple2Equal(tuple->data()[i-1], tuple->data()[i])) {
                    B->data()[i] = B->data()[i-1];
                }
            else {
                B->data()[i] = i;
            }
        }

        indexOffset = lastPossibleIndex + 1;
        lastTuple = tuple->data()[lastPossibleIndex];

        if (worldSize > 1) {
            MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, 1, root, MPI_COMM_WORLD);
            MPI_Send(&lastTuple, 1, MPI_Tuple2, 1, root, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(&indexOffset, 1,         MPI_LONG_LONG_INT, rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&lastTuple, 1,           MPI_Tuple2,        rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int64 size = tuple->size();
        int64 lastPossibleIndex = size-1;
        B->resize(size);
        B->data()[0] = tuple2Equal(lastTuple, tuple->data()[0]) ? indexOffset-1 : indexOffset;

        for (int64 i = 1; i < size; i++) {
            if (tuple2Equal(tuple->data()[i-1], tuple->data()[i])) {
                    B->data()[i] = B->data()[i-1];
                }
            else {
                B->data()[i] = i + indexOffset;
            }
        }

        indexOffset += lastPossibleIndex + 1;
        lastTuple = tuple->data()[lastPossibleIndex];

        if (rank < worldSize-1) {
            MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, rank+1, rank, MPI_COMM_WORLD);
            MPI_Send(&lastTuple, 1, MPI_Tuple2, rank+1, rank, MPI_COMM_WORLD);
        }
    }

}
