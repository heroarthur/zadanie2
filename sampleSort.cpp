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

#include "sortsOpenMP.cpp"

#include<bits/stdc++.h>


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

bool tuple3Equal(Tuple3 t1, Tuple3 t2) {
    return t1.B == t2.B && t1.B2 == t2.B2;
}

bool tuple3Greater(Tuple3 t1, Tuple3 t2) {
    return t1.B > t2.B || (t1.B == t2.B && t1.B2 > t2.B2);
}

bool tuple3Smaller(Tuple3 t1, Tuple3 t2) {
    return !tuple3Equal(t1, t2) && !tuple3Greater(t1, t2);
}


// int binarySearchTuple3(Tuple3* arr, Tuple3 tuple, int64 l, int64 r)
// {
//     if (r >= l) {
//         int mid = l + (r - l) / 2;

//         if (tuple3Equal(arr[mid], tuple) || (tuple3Smaller(tuple, arr[mid]) && tuple3Greater(tuple, arr[mid-1])))
//             return mid;

//         if (tuple3Greater(arr[mid], tuple))
//             return binarySearch(arr, l, mid - 1, x);
 
//         return binarySearch(arr, mid + 1, r, x);
//     }
//     return -1;
// }



void sample_sort_MPI_tuple3(vector<Tuple3>* A, 
                            vector<Tuple3>* sample,
                            vector<Tuple3>* rootSampleRecv,
                            vector<Tuple3>* broadcastSample,
                            int64* sizeOnNode, 
                            int rank, 
                            int worldSize) {
    const int root = 0;

    local_sort_openMP_tuple3(A, *sizeOnNode);
    cout<<"local sorted "<<worldSize<<" "<<rank<<endl;



    int p2 = worldSize * worldSize;


    int64 step = *sizeOnNode / (worldSize + 1);
    
    int sendNumber = worldSize;
    for (int i = 0; i < worldSize; i++) {
        sample->data()[i] = A->data()[(i + 1) * step];
        cout<<"sample " << sample->data()[i].B <<", " << sample->data()[i].B2 << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(sample->data(), sendNumber, MPI_Tuple3, (void*)rootSampleRecv->data(), sendNumber, MPI_Tuple3, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // // cout << "gather" << endl;s

    if (rank == root) {
        local_sort_openMP_tuple3(rootSampleRecv, p2);
        for (int i = 0; i < worldSize-1; i++) {
            broadcastSample->data()[i] = rootSampleRecv->data()[(i+1) * worldSize];
        }
    }

    

    MPI_Bcast((void*)broadcastSample->data(), worldSize-1, MPI_Tuple3, root, MPI_COMM_WORLD);



    // for (int i = 0; i < worldSize; i++) {
        // MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums, num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // }

    // wyslij

}








// void sample_sort_MPI_Tuple3() {

// }