#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>

#define int64 long long int
#define K 33
#define ENTER cout<<endl<<endl


using namespace std;




struct Tuple2 {
    char B[K];
    int64 i;
};

// struct Tuple3 {
//     int64 B;
//     int64 B2;
//     int64 i;
// };

typedef struct mpi_tuple3 {
    int64 B;
    int64 B2;
    int64 i;
} Tuple3;


MPI_Datatype MPI_Tuple3;

int get_block_start(int blockId, int blocksNumber, int dataSize) {
	return (dataSize / blocksNumber) * blockId;
}


// void local_sort_openMP_tuple2(Tuple2* A, int64 size, int k) {

//     auto cmp = [](Tuple2 a, Tuple2 b) { 
//         return strcmp(a.B, b.B) >= 0; 
//     };

// 	int blocksNumber = 0;
// 	#pragma omp parallel
// 	{
// 		blocksNumber = omp_get_num_threads();
// 		int lastElemensSize = size % blocksNumber;
// 		int blockId = omp_get_thread_num();
// 		int blockStart = get_block_start(blockId, blocksNumber, size);
// 		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
// 		std::sort(A + blockStart, A + blockEnd, cmp);
// 	}

// 	int merges = blocksNumber / 2;
// 	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
// 	{
// 		int mergesInStep = (blocksNumber / (2 * mergeStep));

// 		#pragma omp parallel for
// 		for (int i = 0; i < mergesInStep; i++) {
// 			int64 halfMergeLen = (size / blocksNumber) * mergeStep;
// 			int64 mergeStart = i * 2 * halfMergeLen;
// 			int64 mergeMid = mergeStart + halfMergeLen;
// 			int64 mergeEnd = i == mergesInStep-1 ? size : mergeStart + 2 * halfMergeLen;
// 			mergeEnd = min(mergeEnd, size);
// 			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, cmp);
// 		}
// 	}
// }





void local_sort_openMP_tuple3(vector<Tuple3>& A, int64 size) {
    
	struct cmp_tuple3 {
		bool operator ()(Tuple3 const& a, Tuple3 const& b) const {
			return a.B < b.B || (a.B == b.B && a.B2 < b.B2);
		}
	};

    // posortowac wszystkie czesci lokalnie
	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = size % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, size);
		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A + blockStart, A + blockEnd, cmp_tuple3());
	}

	int merges = blocksNumber / 2;
	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
	{
		int mergesInStep = (blocksNumber / (2 * mergeStep));

		#pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
			int64 halfMergeLen = (size / blocksNumber) * mergeStep;
			int64 mergeStart = i * 2 * halfMergeLen;
			int64 mergeMid = mergeStart + halfMergeLen;
			int64 mergeEnd = i == mergesInStep-1 ? size : mergeStart + 2 * halfMergeLen;
			mergeEnd = min(mergeEnd, size);
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, cmp_tuple3());
		}
	}
}


