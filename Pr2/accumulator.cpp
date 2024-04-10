#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <omp.h>
#include <cstdlib>
#include <ctime>

using namespace std;


//ArowVal, ArowInd, ArowPtr, BrowVal, BrowInd, BrowPtr, size
double* multiplyCSRMatrix(double* ArowVal, double* ArowInd, double* ArowPtr, double* BrowVal, double* BrowInd, double* BrowPtr, int ArowPtrSize, int BrowPtrSize, int A_rows, int B_cols, int numThreads) {
    double* result = new double[A_rows * B_cols];
    for(int i = 0; i < A_rows * B_cols; i++) {
        result[i] = 0;
    }
    int maxThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    double start;
    double end;
    double temp;

    
    int count = 0;
    double temp1, temp2;

    start = omp_get_wtime();
    
    #pragma omp parallel for shared(result)
    for(int i = 0; i < ArowPtrSize - 1; i++) {
        for(int j = 0; j < BrowPtrSize - 1; j++) {
            for(int k = ArowPtr[i]; k < ArowPtr[i + 1]; k++) {
                for(int m = BrowPtr[j]; m < BrowPtr[j + 1]; m++) {
                    temp1 = ArowVal[k];
                    temp2 = BrowVal[m];
                    if(temp1 == temp2){
                        result[count] += temp1 * temp2;
                    } 
                }
            }
            count++;
        }
    }



    end = omp_get_wtime();
    cout << "time par: " << end - start << ", numThreads: " << numThreads << endl;
    //cout << /*"time par: " <<*/ end - start << endl; //<< ", numThreads: " << numThreads << ", total: " << total * numThreads << endl;
    

    return result;
}

double* multiplyNormalMatrix(double* A, double* B, double* buffer, int A_rows, int A_cols, int B_cols, int numThreads) {
    int maxThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    double* result = new double[A_rows * B_cols];
    for(int i = 0; i < A_rows * B_cols; i++) {
        result[i] = 0;
    }

    double start;
    double end;
    start = omp_get_wtime();
    long int total = 0;
    #pragma omp parallel for shared(buffer, result, total)
    for(int i = 0; i < A_rows; ++i) {
        for(int k = 0; k < A_cols; ++k) {
            for(int j = 0; j < B_cols; ++j) {
                buffer[k * B_cols + j] = A[i * A_cols + k] * B[k * B_cols + j];
            }
            for(int j = 0; j < B_cols; ++j) {
                result[i * B_cols + j] += buffer[k * B_cols + j];
                total++;
            }
        }
    }
    end = omp_get_wtime();
    cout << "time par: " << end - start << ", numThreads: " << numThreads << ", total: " << endl;
    //cout << /*"time par: " <<*/ end - start << endl; //<< ", numThreads: " << numThreads << ", total: " << total * numThreads << endl;
    return result;
}



void printMatrix(const double* matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << matrix[i * cols + j] << " ";
        }
        cout << endl;
    }
}


int main(int argc, char* argv[]) {
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt(seed); 
    int min = 0;
    int max = 20;
    uniform_int_distribution<int> dist(min, max);
    int randomNum = dist(mt);

    int nrows = atoi(argv[1]);
    cout << "Input = " << nrows <<endl;
    nrows = atoi(argv[2]);
    cout << "Input = " << nrows <<endl;
    nrows = atoi(argv[3]);
    cout << "Input = " << nrows <<endl;


    double fillPercent = stod(argv[4]);;
    cout << fillPercent << endl;

    srand(static_cast<unsigned int>(time(nullptr)));
    double random_number = static_cast<double>(rand()) / RAND_MAX;
    cout << "random: " << random_number << endl;

    double start;
    double end;

    //define matrices
    //vector<double> A;
    //vector<double> B;
    //vector<double> C;
    
    int A_rows = atoi(argv[1]);
    int A_cols = atoi(argv[2]); // = B_rows
    int B_cols = atoi(argv[3]);;

    int n = A_rows;
    int p = A_cols;
    int m = B_cols;

    //int size = 50;
    double* A = new double[A_rows * A_cols];
    double* B = new double[A_cols * B_cols];
    double* C = new double[A_rows * B_cols];
    
    

    for(int i = 0; i < A_rows * B_cols; i++) {
        C[i] = 0;
    }



    //sparse values
    vector<double> ArowValVec;
    vector<double> ArowIndVec;
    vector<double> ArowPtrVec;
    int ArowPtrSize;
    ArowPtrVec.push_back(0);

    vector<double> BrowValVec;
    vector<double> BrowIndVec;
    vector<double> BrowPtrVec;
    BrowPtrVec.push_back(0);
    int BrowPtrSize;


    int count = 0;
    int row_ptr = 0;
    double temp;
    double sum; 
    int numThreads = 16;
    omp_set_num_threads(numThreads);

    //generate matrices
    //#pragma omp parallel for shared(A) private(i, j, randomNum);
    for(int i = 0; i < A_rows; i++) {
        for(int j = 0; j < A_cols; j++) {
            if(random_number < fillPercent) {
                A[i * A_cols + j] = static_cast<double>(rand()) / RAND_MAX;
            } else {
                A[i * A_cols + j] = 0;
            }
            random_number = static_cast<double>(rand()) / RAND_MAX; 
        }
    }
    //#pragma omp parallel for shared(B) private(i, j, randomNum);
    for(int i = 0; i < A_cols; i++) {
        for(int j = 0; j < B_cols; j++) {
            if(random_number < fillPercent) {
                B[i * B_cols + j] = static_cast<double>(rand()) / RAND_MAX;
            } else {
                B[i * B_cols + j] = 0;
            }
            random_number = static_cast<double>(rand()) / RAND_MAX; 
        }
    }
    //transpose
    for(int i = 0; i < A_cols; ++i) {
        for (int j = i + 1; j < B_cols; ++j) {
           swap(B[i * B_cols + j], B[j * A_cols + i]);
        }
    }
    

    //print matrices
    /*
    cout << "A:" << endl;
    for(int i = 0; i < A_rows; i++) {
        for(int j = 0; j < A_cols; j++) {
            cout << A[i * A_cols + j] << ", ";
        }
        cout << endl;
    }
    cout << endl << "B:" << endl;
    for(int i = 0; i < A_cols; i++) {
        for(int j = 0; j < B_cols; j++) {
            cout << B[i * B_cols + j] << ", ";
        }
        cout << endl;
    }*/
    
    /*
    cout << endl << "B trans:" << endl;
    for(int i = 0; i < A_cols; i++) {
        for(int j = 0; j < B_cols; j++) {
            cout << B[i * B_cols + j] << ", ";
        }
        cout << endl;
    }*/

    //update sparce row values
    for(int i = 0; i < A_rows; i++) {
        for(int j = 0; j < A_cols; j++) {
            temp = A[i * A_cols + j];
            if(temp != 0) {
                ArowValVec.push_back(temp);
                ArowIndVec.push_back(j);
            }
        }
        ArowPtrVec.push_back(ArowValVec.size());
    }
    for(int i = 0; i < A_cols; i++) {
        for(int j = 0; j < B_cols; j++) {
            temp = B[i * B_cols + j];
            if(temp != 0) {
                BrowValVec.push_back(temp);
                BrowIndVec.push_back(j);
            }
        }
        BrowPtrVec.push_back(BrowValVec.size());
    }
    //create non-vector arrays for sparse values
    double* ArowVal = new double[ArowValVec.size()];
    double* ArowInd = new double[ArowIndVec.size()];
    double* ArowPtr = new double[ArowPtrVec.size()];

    double* BrowVal = new double[BrowValVec.size()];
    double* BrowInd = new double[BrowIndVec.size()];
    double* BrowPtr = new double[BrowPtrVec.size()];

    for(int i = 0; i < ArowValVec.size(); i++) {
        ArowVal[i] = ArowValVec[i];
        ArowInd[i] = ArowIndVec[i];
    }
    for(int i = 0; i < ArowPtrVec.size(); i++) {
        ArowPtr[i] = ArowPtrVec[i];
    }
    for(int i = 0; i < BrowValVec.size(); i++) {
        BrowVal[i] = BrowValVec[i];
        BrowInd[i] = BrowIndVec[i];
    }
    for(int i = 0; i < BrowPtrVec.size(); i++) {
        BrowPtr[i] = BrowPtrVec[i];
    }
    ArowPtrSize = ArowPtrVec.size();
    BrowPtrSize = BrowPtrVec.size();

    //prints rowval, worind, and rowptr
    /*
    cout << "A row val: ";
    for(int i = 0; i < ArowVal.size(); i++) {
        cout << ArowVal[i] << ", ";
    }
    cout << endl << "A row ind: ";
    for(int i = 0; i < ArowInd.size(); i++) {
        cout << ArowInd[i] << ", ";
    }
    cout << endl << "A row ptr: ";
    for(int i = 0; i < ArowPtr.size(); i++) {
        cout << ArowPtr[i] << ", ";
    }
    cout << endl;

    cout << "B row val: ";
    for(int i = 0; i < BrowVal.size(); i++) {
        cout << BrowVal[i] << ", ";
    }
    cout << endl << "B row ind: ";
    for(int i = 0; i < BrowInd.size(); i++) {
        cout << BrowInd[i] << ", ";
    }
    cout << endl << "B row ptr: ";
    for(int i = 0; i < BrowPtr.size(); i++) {
        cout << BrowPtr[i] << ", ";
    }
    cout << endl;
    */

    //multiplication
    double* buffer = new double[A_cols * B_cols];

    
    /*
    for(int i = 0; i < A_rows; i++) {
        for(int j = 0; j < B_cols; j++) {
            cout << C[i * B_cols + j] << ", ";
        }
        cout << endl;
    }*/
    
    
    // Initialize the buffer matrix with zeros
    int i, j, k;
    // Perform matrix multiplication
    cout << "Normal multiplication method in parallel: " << endl;
    for(int threads = 28; threads < 29; threads++) {
        if(threads % 4 == 0 || threads == 1 || threads == 2) {
            omp_set_num_threads(threads);
            //start = omp_get_wtime();
            C = multiplyNormalMatrix(A, B, buffer, A_rows, A_cols, B_cols, threads);
            //end = omp_get_wtime();
            //cout << "time par: " << end - start << ", numThreads: " << threads << endl;
        }
    }
    //cout << endl << "Result:" << endl;
    //printMatrix(C, A_rows, B_cols);

    /*for(int i = 0; i < A_rows * B_cols; i++) {
        C[i] = 0;
    }*/
    cout << "Accumulator method in parallel: " << endl;
    for(int threads = 28; threads < 29; threads++) {
        if(threads % 4 == 0 || threads == 1 || threads == 2) {
            //omp_set_num_threads(threads);
            /*start = omp_get_wtime();
            #pragma omp parallel for shared(C) //private(i, j, k, kA, kB, colA, colB)
            for(int i = 0; i < size; ++i) {
                for(int j = 0; j < size; ++j) {
                    for(int kA = ArowPtr[i]; kA < ArowPtr[i + 1]; ++kA) {
                        int colA = ArowInd[kA];
                        for(int kB = BrowPtr[colA]; kB < BrowPtr[colA + 1]; ++kB) {
                            int colB = BrowInd[kB];
                            C[i * size + colB] += ArowVal[kA] * BrowVal[kB];
                        }
                    }
                }
            }*/
            //start = omp_get_wtime();
            omp_set_num_threads(threads);
            C = multiplyCSRMatrix(ArowVal, ArowInd, ArowPtr, BrowVal, BrowInd, BrowPtr, ArowPtrSize, BrowPtrSize, A_rows, B_cols, threads);
            //end = omp_get_wtime();
            //cout << "time par: " << end - start << ", numThreads: " << threads << endl;
        }
    }
    //cout << endl << "Result:" << endl;
    //printMatrix(C, A_rows, B_cols);

    return 0;
}



