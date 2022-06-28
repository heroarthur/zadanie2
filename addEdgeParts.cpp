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


#include<bits/stdc++.h>

#define root 0


using namespace std;


// int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
void addEdgeParts(vector<char>* __restrict__ nodeCharArray, 
                  int rank, 
                  int worldSize) {
    
    vector<char> receivedEdgePart, sendEdgePart;



    if (rank == worldSize-1) {
        sendEdgePart.resize(minInt64(K-1, nodeCharArray->size()));
		copy(nodeCharArray->begin(), nodeCharArray->begin() + minInt64(K-1, nodeCharArray->size()), sendEdgePart.data());

        if (worldSize > 1) {
            MPI_Send(sendEdgePart.data(), sendEdgePart.size(), MPI_CHAR, rank-1, rank, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Status status;
        int receivedSize;
        receivedEdgePart.resize(K-1);
        MPI_Recv(receivedEdgePart.data(), K-1, MPI_CHAR, rank+1, rank+1, MPI_COMM_WORLD, &status); 
        MPI_Get_count(&status, MPI_CHAR, &receivedSize);

        nodeCharArray->insert(nodeCharArray->end(), receivedEdgePart.begin(), receivedEdgePart.begin() + receivedSize);
        sendEdgePart.resize(minInt64(K-1, nodeCharArray->size()));
        copy(nodeCharArray->begin(), nodeCharArray->begin() + minInt64(K-1, nodeCharArray->size()), sendEdgePart.data());

        if (rank > 0) {
            MPI_Send(sendEdgePart.data(), sendEdgePart.size(), MPI_CHAR, rank-1, rank, MPI_COMM_WORLD);
        }
    }
    
}