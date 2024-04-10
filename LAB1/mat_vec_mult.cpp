#include <iostream>
#include <chrono>
#include "Matrix.h"


int main(int argc, char **argv){
    Matrix mat;
    Vector vec;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> res = mat * vec;
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("Serial Execution Time (microseconds): %ld\n", duration.count());
    return 0;
}
