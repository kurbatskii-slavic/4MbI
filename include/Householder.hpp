#include "matrix.hpp"
#include <vector>

void
QRHouseholder(std::vector<HouseholderMatrix> &H, Matrix &A) // QR decomposition using Householder method
{
    HouseholderMatrix T;
    std::vector<double> ref(A[0].size());
    for (size_t i = 0; i < A.cols - 1; i++) {
        get_reflection(ref, A[i], i);
        make_reflection(T, A[i], ref, i);
        H.push_back(T);
        for (size_t j = i; j < A.cols; j++) {
            matvec(T, A[j]);
        }
    }        
}

std::ostream 
&operator<<(std::ostream &os, std::pair<std::vector<HouseholderMatrix>, Matrix> &result) // cout overloading for result pair
{
    for (auto T: result.first) { // display H
        os << T << std::endl;
    }
    os << std::endl << result.second << std::endl; // display R
    return os;
}
