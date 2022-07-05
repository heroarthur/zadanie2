#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <fstream> 

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include "SA.cpp"

using namespace std;



int main(int argc, char** argv) {

	srand(1);

    // Read from the text file
    string fileName(argv[1]);
    ifstream MyReadFile(fileName);
    string genomeArray;
    getline(MyReadFile, genomeArray);
    MyReadFile.close();

    vector<char> inputGenome;
    inputGenome.insert(inputGenome.begin(), genomeArray.c_str(), genomeArray.c_str() + genomeArray.size());
    inputGenome.push_back('$');

    int64 dataSize = inputGenome.size();
	
	vector<int64> B_1;
	vector<int64> B_2;

	vector<Tuple2> tuple2_Arr; tuple2_Arr.resize(dataSize); 
	vector<Tuple2> tuple2_second; 

	vector<Tuple3> tuple3; 
	vector<Tuple3> tuple3_second; 

	vector<int64> SA, SA_second;

	bool allSingletones;

    for (int i = 0; i < dataSize; i++) {
        memset(tuple2_Arr[i].B, 0, K);
        copy(inputGenome.begin() + i, inputGenome.begin() + i + minInt64(k, dataSize), tuple2_Arr[i].B);
        tuple2_Arr[i].i = i;
    }



    
    std::sort(tuple2_Arr.begin(), tuple2_Arr.end(), cmp_tuple2());

    B_1.resize(dataSize);
    B_1.data()[0] = 0;

    B_2.resize(dataSize);
    tuple3.resize(B_1.size());

    for (int64 i = 1; i < dataSize; i++) {
        if (tuple2Equal(tuple2_Arr.data()[i-1], tuple2_Arr.data()[i])) {
                B_1.data()[i] = B_1.data()[i-1];
            }
        else {
            B_1.data()[i] = i;
        }
        // cout<<tuple2_Arr.data()[i].B<<" ";
    }
    // cout<<endl;

    SA.resize(tuple2_Arr.size());
    SA_second.resize(tuple2_Arr.size());


    for (int i = 0; i < tuple2_Arr.size(); i++) {
        SA.data()[i] = tuple2_Arr.data()[i].i;
    }


    bool done = false;
	for (int64 h = k; true; h*=2) {

		// print_vector(&B_1);
		// print_vector(&SA);
        // cout<<endl;

        for (int i = 0; i < dataSize; i++) {
            B_2[SA[i]] = B_1[i];
            SA_second[SA[i]] = SA[i];
        }
        for (int i = 0; i < dataSize; i++) {
            B_1[i] = B_2[i];
            // SA[i] = SA_second[i];
        }

        if (done) {
            break;
        }

        for (int i = 0; i < dataSize; i++) {
            SA[i] = SA_second[i];
        }

        std::fill(B_2.begin(), B_2.end(), 0);

        for (int64 i = 0; i + h < dataSize; i++) {
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

        // for (int i = 0; i < dataSize; i++) {
        //     cout<<"("<<tuple3.data()[i].B<<" "<<tuple3.data()[i].i<<<<") ";
        // }
        // cout<<endl<<endl;

        for (int64 i = 1; i < dataSize; i++) {
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


