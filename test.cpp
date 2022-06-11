#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>
// #include <mpi.h>


#define int64 long long int
#define root 0

using namespace std;





// void sample_sort_MPI()

// struct Tuple2 {
//     char B[K];
//     int64 i;
// };

// struct Tuple3 {
//     int64 B;
//     int64 B2;
//     int64 i;
// };



int main(int argc, char** argv) {

    int v;
    // #pragma omp parallel for
    for (int i = 0; i < 10000000; i++) {
        v++;
    }

    cout<<v<<endl;
	return 0;
}


