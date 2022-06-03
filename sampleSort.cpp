#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>

#include <omp.h>
#include <mpi.h>

#include "sortsOpenMP.cpp"


using namespace std;

// struct Tuple2 {
//     char B[K];
//     int64 i;
// };

// struct Tuple3 {
//     int64 B;
//     int64 B2;
//     int64 i;
// };








void sample_sort_MPI_tuple3(Tuple3* A, 
                            Tuple3* sample,
                            Tuple3* rootSampleRecv,
                            Tuple3* broadcastSample,
                            int64 sizeOnNode, 
                            int rank, 
                            int worldSize) {
    const int root = 0;

    local_sort_openMP_tuple3(A, sizeOnNode);
    cout<<"local sorted "<<worldSize<<" "<<rank<<endl;

    // if (rank == 3) {
	// 	for (int i = 0; i < sizeOnNode; i++) {
	// 		cout<<"("<<A[i].B<<","<<A[i].B2<<") ";
	// 	}
	// }
	// ENTER;

    int p2 = worldSize * worldSize;


    int64 step = sizeOnNode / (worldSize + 1);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        sample[i] = A[(i + 1) * step];
        cout<<"sample " << sample[i].B <<", " << sample[i].B2 << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&sample, sendNumber, mpi_tuple3, rootSampleRecv, sendNumber, mpi_tuple3, root, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    // cout << "gather" << endl;s

    if (rank == root) {


        if (rank == 0) {
        	for (int i = 0; i < p2; i++) {
        		cout<<"("<<rootSampleRecv[i].B<<","<<rootSampleRecv[i].B2<<") ";
        	}
        }
        // local_sort_openMP_tuple3(rootSampleRecv, p2);
        ENTER;

        // for (int i = 0; i < worldSize-1; i++) {
        //     broadcastSample[i] = rootSampleRecv[(i+1) * worldSize];
        //     cout<<"("<<broadcastSample[i].B<<", "<<broadcastSample[i].B2<<")"<<endl;
        // }
    }

    

    // MPI_Bcast(broadcastSample, worldSize-1, MPI_LONG_LONG_INT, root, MPI_COMM_WORLD);

    // if (rank == root) {
    //     for (int i = 0; i < worldSize-1; i++) {
    //         cout<<"<"<<broadcastSample[i].B<<","<<broadcastSample[i].B2<<"> ";
    //     }
    //     ENTER;
    // }

}








// void sample_sort_MPI_tuple3() {

// }