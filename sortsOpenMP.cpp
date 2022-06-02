#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <omp.h>

#define int64 long long int
#define K 33


using namespace std;




struct Tuple2 {
    char B[K];
    int64 i;
};

struct Tuple3 {
    int64 B;
    int64 B2;
    int64 i;
};

int get_block_start(int blockId, int blocksNumber, int dataSize) {
	return (dataSize / blocksNumber) * blockId;
}

void local_sort_openMP_tuple2(Tuple2* A, int size, int k) {

    auto cmp = [](Tuple2 a, Tuple2 b) { 
        return strcmp(a.B, b.B) >= 0; 
    };

	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = size % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, size);
		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A + blockStart, A + blockEnd, cmp);
	}

	int merges = blocksNumber / 2;
	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
	{
		int mergesInStep = (blocksNumber / (2 * mergeStep));

		#pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
			int halfMergeLen = (size / blocksNumber) * mergeStep;
			int mergeStart = i * 2 * halfMergeLen;
			int mergeMid = mergeStart + halfMergeLen;
			int mergeEnd = i == mergesInStep-1 ? size : mergeStart + 2 * halfMergeLen;
			mergeEnd = min(mergeEnd, size);
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, cmp);
		}
	}
}



void local_sort_openMP_tuple3(Tuple3* A, int size) {

    auto cmp = [](Tuple3 a, Tuple3 b) {
        if (a.B > b.B) {
            return true;
        }
        else if (a.B2 > b.B2) {
            return true;
        }
        return false;
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
		std::sort(A + blockStart, A + blockEnd, cmp);
	}

	int merges = blocksNumber / 2;
	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
	{
		int mergesInStep = (blocksNumber / (2 * mergeStep));

		#pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
			int halfMergeLen = (size / blocksNumber) * mergeStep;
			int mergeStart = i * 2 * halfMergeLen;
			int mergeMid = mergeStart + halfMergeLen;
			int mergeEnd = i == mergesInStep-1 ? size : mergeStart + 2 * halfMergeLen;
			mergeEnd = min(mergeEnd, size);
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, cmp);
		}
	}
}


