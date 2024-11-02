#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <chrono>
#include <random>


#include "matrix.hpp"
#include "Householder.hpp"
#include "Gramm-Shmidt.hpp"
#include "generate_rand.hpp"


enum
{
    SIZE = 100 // from the variant
};

int
main()
{
    Matrix A(SIZE, SIZE); 
    std::cin >> A; // read matrix from file
    std::vector<HouseholderMatrix> H; // vector for Householder matrices
    Matrix A_f = A;
    std::vector<double> x = generate_random(-1, 1, SIZE); // generate random vector
    std::vector<double> f = A * x; // right side of Ax = f system
    #ifdef TIME // time measure
        const auto start{std::chrono::steady_clock::now()}; // start time
    #endif
    QRHouseholder(H, A_f); // solve QR decomposition with Householder method
    #ifdef TIME
        const auto end{std::chrono::steady_clock::now()}; // end time
        const std::chrono::duration<double> elapsed_seconds{end - start}; // duration
        std::cout << "Time: " << elapsed_seconds.count() * 1000 << "ms" << std::endl; // display duration
    #endif
    for (size_t j = 0; j < A_f.cols; j++) {
        for (size_t i = j + 1; i < A_f.rows; i++) {
            A_f(i, j) = 0; // place 0 for better precision
        }
    }
    std::cout << A_f << '\n';
    std::pair<std::vector<HouseholderMatrix>, Matrix> result = std::make_pair(H, A_f); // store the result
    Matrix R = A_f;
    for (int i = H.size() - 1; i >= 0; i--) {
        for (size_t j = 0; j < R.cols; j++) {
            matvec(H[i], R[j]);
        }
    }
    std::cout << "Householder: ||A - QR|| = " << matrix_maximum_norm(A - R) << std::endl; 
    std::vector<double> f_n = f;
    for (auto &T: H) matvec(T, f_n); // transform Ax = f -> Rx = Q^T f
    std::vector<double> x_f(SIZE);
    x_f = solve_triangular_system(A_f, f_n);
    double f_error = maximum_norm(f - A * x_f); // calculate f error 
    double x_error = maximum_norm(x - x_f); // calculate x error
    std::cout << "||x - x_f|| = " <<  x_error << std::endl;
    std::cout << "||f - Ax_f|| = " << f_error << std::endl;
 
    return 0;
}