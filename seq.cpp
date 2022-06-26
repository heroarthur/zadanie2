#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include "SA.cpp"

using namespace std;



int main(int argc, char** argv) {

	srand(1);

	int64 allDataSize = 100;
	int64 singleNodeDataSize = allDataSize;
	
	vector<int64> B_1;
	vector<int64> B_2;

	vector<Tuple2> tuple2_Arr; tuple2_Arr.resize(singleNodeDataSize); 
	vector<Tuple2> tuple2_second; 

	vector<Tuple3> tuple3; 
	vector<Tuple3> tuple3_second; 





	vector<int64> SA, SA_second;

	bool allSingletones;

    for (int i = 0; i < singleNodeDataSize; i++) {
        fillCharArray(tuple2_Arr[i].B);
        tuple2_Arr[i].i = i;
    }
    tuple2_Arr[singleNodeDataSize-1].B[charArrayLen-1] = '$';
    
    std::sort(tuple2_Arr.begin(), tuple2_Arr.end(), cmp_tuple2());


    B_1.resize(singleNodeDataSize);
    B_1.data()[0] = 0;

    B_2.resize(singleNodeDataSize);
    tuple3.resize(B_1.size());

    for (int64 i = 1; i < singleNodeDataSize; i++) {
        if (tuple2Equal(tuple2_Arr.data()[i-1], tuple2_Arr.data()[i])) {
                B_1.data()[i] = B_1.data()[i-1];
            }
        else {
            B_1.data()[i] = i;
        }
    }

    bool done = false;
	for (int64 h = k; true; h*=2) {

        SA.resize(tuple2_Arr.size());
        for (int i = 0; i < tuple2_Arr.size(); i++) {
            SA.data()[i] = tuple2_Arr.data()[i].i;
        }


        std::fill(B_2.begin(), B_2.end(), -1);

        if (done) {
            break;
        }

        for (int i = 0; i + h < singleNodeDataSize; i++) {
            B_2[i] = B_1[i+h];
        }

        for (int i = 0; i < tuple3.size(); i++) {
            tuple3.data()[i].B = B_1.data()[i];
            tuple3.data()[i].B2 = B_2.data()[i];
            tuple3.data()[i].i = SA.data()[i];
        }

		std::sort(tuple3.begin(), tuple3.end(), cmp_tuple3());

        done = true;
		B_1.data()[0] = 0;
        SA.data()[0] = tuple3.data()[0].i;
        for (int64 i = 1; i < singleNodeDataSize; i++) {
            if (tuple3Equal(tuple3.data()[i-1], tuple3.data()[i])) {
                    done = false;
                    B_1.data()[i] = B_1.data()[i-1];
                }
            else {
                B_1.data()[i] = i;
            }
            SA.data()[i] = tuple3.data()[i].i;
        }
	}

    print_vector(&SA);

	return 0;
}


