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

#define wyslijRaz 5

using namespace std;



void rebucketing_2h_group_rank(vector<Tuple3>* tuple, 
                               vector<int64>* B, 
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
        B->data()[0] = 0;
        for (int64 i = 1; i < size; i++) {
            if (tuple3Equal(tuple->data()[i-1], tuple->data()[i])) {
                    localAllSingletones = false;
                    B->data()[i] = B->data()[i-1];
                }
            else {
                B->data()[i] = i;
            }
        }

        allSingletonesCheck = localAllSingletones;
        indexOffset = lastPossibleIndex;
        lastTuple = tuple->data()[lastPossibleIndex];

        MPI_Send(&allSingletonesCheck, 1, MPI_C_BOOL, 1, root, MPI_COMM_WORLD);
        MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, 1, root, MPI_COMM_WORLD);
        MPI_Send(&lastTuple, 1, MPI_Tuple3, 1, root, MPI_COMM_WORLD);
    }
    else {
        MPI_Recv(&allSingletonesCheck, 1, MPI_C_BOOL,        rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&indexOffset, 1,         MPI_LONG_LONG_INT, rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&lastTuple, 1,           MPI_Tuple3,        rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int64 size = tuple->size();
        int64 lastPossibleIndex = size-1;
        B->resize(size);
        B->data()[0] = 0;
        for (int64 i = 1; i < size; i++) {
            if (tuple3Equal(tuple->data()[i-1], tuple->data()[i])) {
                    allSingletonesCheck = false;
                    B->data()[i] = B->data()[i-1];
                }
            else {
                B->data()[i] = i + indexOffset;
            }
        }

        indexOffset += lastPossibleIndex;
        lastTuple = tuple->data()[lastPossibleIndex];

        if (rank < worldSize-1) {
            MPI_Send(&allSingletonesCheck, 1, MPI_C_BOOL, rank+1, rank, MPI_COMM_WORLD);
            MPI_Send(&indexOffset, 1, MPI_LONG_LONG_INT, rank+1, rank, MPI_COMM_WORLD);
            MPI_Send(&lastTuple, 1, MPI_Tuple3, rank+1, rank, MPI_COMM_WORLD);
        }
    }

    MPI_Bcast(&allSingletonesCheck, 1, MPI_C_BOOL, worldSize-1, MPI_COMM_WORLD);
}


