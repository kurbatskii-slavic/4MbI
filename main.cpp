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
#include "Gram-Shmidt.hpp"
#include "generate_rand.hpp"


enum
{
    SCALE = 1000,
    BASE = 10
};

int
main(int argc, char** argv)
{
    const size_t SIZE = std::strtol(argv[1], nullptr, BASE);
    Matrix A(SIZE, SIZE); 
    std::cin >> A; // read matrix from file
    Matrix R_h = A; // save A for future
    std::vector<HouseholderMatrix> H; // vector for Householder matrices
    Matrix R(SIZE, SIZE), Q(SIZE, SIZE); // matrices for Gramm-Shmidt algorithm
    std::vector<double> x = generate_random(-1, 1, SIZE); // generate random vector
    std::vector<double> f = A * x; // right side of Ax = f system
    #ifdef TIME // time measure
        const auto start{std::chrono::steady_clock::now()}; // start time
    #endif
    QRHouseholder(H, R_h); // solve QR decomposition with Householder method
    #ifdef TIME
        const auto end{std::chrono::steady_clock::now()}; // end time
        const std::chrono::duration<double> elapsed_seconds{end - start}; // duration
        std::cout << "Time (Householder): " << elapsed_seconds.count() * SCALE << "ms" << std::endl; // display duration
    #endif
    for (size_t j = 0; j < R_h.cols; j++) {
        for (size_t i = j + 1; i < R_h.rows; i++) {
            R_h(i, j) = 0; // place 0 for better precision
        }
    }
    Matrix R_f = R_h;
    for (int i = H.size() - 1; i >= 0; i--) {
        for (size_t j = 0; j < R_f.cols; j++) {
            matvec(H[i], R_f[j]); // calculate QR
        }
    }
    std::cout << "||A - QR|| = " << matrix_maximum_norm(A - R_f) << std::endl << std::endl; 
    
    #ifdef TIME // time measure
        const auto start_2{std::chrono::steady_clock::now()}; // start time
    #endif
    QRGram_Shmidt(Q, R, A);
    #ifdef TIME
        const auto end_2{std::chrono::steady_clock::now()}; // end time
        const std::chrono::duration<double> elapsed_seconds_2{end_2 - start_2}; // duration
        std::cout << "Time: (Gramm-Shmidt): " << elapsed_seconds_2.count() * SCALE << "ms" << std::endl; // display duration
    #endif
    for (size_t j = 0; j < R.cols; j++) {
        for (size_t i = j + 1; i < R.rows; i++) {
            R(i, j) = 0; // place 0 for better precision
        }
    }
    std::cout << "||A - QR|| = " << matrix_maximum_norm(A - Q * R) << std::endl << std::endl; 


    // Householder
    std::vector<double> f_n = f;
    std::vector<double> x_f(SIZE);
    #ifdef TIME // time measure
        const auto start_4{std::chrono::steady_clock::now()}; // start time
    #endif
    for (auto &T: H) matvec(T, f_n); // transform Ax = f -> Rx = Q^T f
    x_f = solve_triangular_system(R_h, f_n);
    #ifdef TIME
        const auto end_4{std::chrono::steady_clock::now()}; // end time
        const std::chrono::duration<double> elapsed_seconds_4{end_4 - start_4}; // duration
        std::cout << "Time: (Householder system): " << elapsed_seconds_4.count() * SCALE << "ms" << std::endl; // display duration
    #endif
    double f_error = maximum_norm(f - A * x_f); // calculate f error 
    double x_error = maximum_norm(x - x_f); // calculate x error
    std::cout << "||x - x_f|| = " <<  x_error << std::endl;
    std::cout << "||f - Ax_f|| = " << f_error << std::endl << std::endl;

    //Gramm-Shmidt
    #ifdef TIME // time measure
        const auto start_3{std::chrono::steady_clock::now()}; // start time
    #endif
    Q.transpose();
    f_n = Q * f;
    x_f = solve_triangular_system(R, f_n);
    #ifdef TIME
        const auto end_3{std::chrono::steady_clock::now()}; // end time
        const std::chrono::duration<double> elapsed_seconds_3{end_3 - start_3}; // duration
        std::cout << "Time: (Gramm-Shmidt system): " << elapsed_seconds_3.count() * SCALE << "ms" << std::endl; // display duration
    #endif
    f_error = maximum_norm(f - A * x_f); // calculate f error 
    x_error = maximum_norm(x - x_f); // calculate x error
    std::cout << "||x - x_f|| = " <<  x_error << std::endl;
    std::cout << "||f - Ax_f|| = " << f_error << std::endl;
    return 0;
}