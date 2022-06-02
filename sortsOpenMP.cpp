#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#define int64 long long int

using namespace std;




struct Tuple2 {
    char* B;
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

void local_sort_openMP(int64* A, int size) {

	// posortowac wszystkie czesci lokalnie
	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = size % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, size);
		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A + blockStart, A + blockEnd, less<int64>());  //tutaj dodac min zeby nie wywalo 
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
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, less<int64>());
		}
	}
}



void local_sort_openMP(int64* A, int size) {

	// posortowac wszystkie czesci lokalnie
	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = size % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, size);
		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A + blockStart, A + blockEnd, less<int64>());  //tutaj dodac min zeby nie wywalo 
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
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, less<int64>());
		}
	}
}




void local_sort_openMP(int64* A, int size) {

	// posortowac wszystkie czesci lokalnie
	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = size % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, size);
		int blockEnd = get_block_start(blockId+1, blocksNumber, size) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A + blockStart, A + blockEnd, less<int64>());  //tutaj dodac min zeby nie wywalo 
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
			inplace_merge(A + mergeStart, A + mergeMid, A + mergeEnd, less<int64>());
		}
	}
}
