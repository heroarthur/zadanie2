#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>

#define K 13
#define ENTER cout<<endl<<endl


using namespace std;


#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif











void local_sort_openMP_tuple3(vector<Tuple3>* A) {
    
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
		int lastElemensSize = A->size() % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, A->size());
		int blockEnd = get_block_start(blockId+1, blocksNumber, A->size()) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A->begin() + blockStart, A->begin() + blockEnd, cmp_tuple3());
	}

	int merges = blocksNumber / 2;
	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
	{
		int mergesInStep = (blocksNumber / (2 * mergeStep));

		#pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
			int64 halfMergeLen = (A->size() / blocksNumber) * mergeStep;
			int64 mergeStart = i * 2 * halfMergeLen;
			int64 mergeMid = mergeStart + halfMergeLen;
			int64 mergeEnd = i == mergesInStep-1 ? A->size() : mergeStart + 2 * halfMergeLen;
			mergeEnd = minInt64(mergeEnd, A->size());
			inplace_merge(A->begin() + mergeStart, A->begin() + mergeMid, A->begin() + mergeEnd, cmp_tuple3());
		}
	}
}




void local_sort_openMP_tuple2(vector<Tuple2>* A) {
    
	struct cmp_tuple2 {
		bool operator ()(Tuple2 const& a, Tuple2 const& b) const {
			return strcmp(a.B, b.B) < 0;
		}
	};

    // posortowac wszystkie czesci lokalnie
	int blocksNumber = 0;
	#pragma omp parallel
	{
		blocksNumber = omp_get_num_threads();
		int lastElemensSize = A->size() % blocksNumber;
		int blockId = omp_get_thread_num();
		int blockStart = get_block_start(blockId, blocksNumber, A->size());
		int blockEnd = get_block_start(blockId+1, blocksNumber, A->size()) + (blockId == blocksNumber-1 ? lastElemensSize : 0);
		std::sort(A->begin() + blockStart, A->begin() + blockEnd, cmp_tuple2());
	}

	int merges = blocksNumber / 2;
	for (int mergeStep = 1; mergeStep < blocksNumber; mergeStep *= 2 )
	{
		int mergesInStep = (blocksNumber / (2 * mergeStep));

		#pragma omp parallel for
		for (int i = 0; i < mergesInStep; i++) {
			int64 halfMergeLen = (A->size() / blocksNumber) * mergeStep;
			int64 mergeStart = i * 2 * halfMergeLen;
			int64 mergeMid = mergeStart + halfMergeLen;
			int64 mergeEnd = i == mergesInStep-1 ? A->size() : mergeStart + 2 * halfMergeLen;
			mergeEnd = minInt64(mergeEnd, A->size());
			inplace_merge(A->begin() + mergeStart, A->begin() + mergeMid, A->begin() + mergeEnd, cmp_tuple2());
		}
	}
}