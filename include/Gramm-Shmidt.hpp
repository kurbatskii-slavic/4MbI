#include "matrix.hpp"
#include <vector>

void
QRGramm_Shmidt(Matrix &Q, Matrix &A) // QR decomposition using Gramm-Shmidt orthogonalization algorithm
{
    for (size_t k = 0; k < Q.cols; k++) {
        std::vector<double> b_k = A[k];
        for (size_t s = 0; s < k; s++) {
            //b_k -= Q[s] * dot_product(A[k], Q[s]);
        }
    }
}