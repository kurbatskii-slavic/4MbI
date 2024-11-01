#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <chrono>

#include "matrix.hpp"
#include "Householder.hpp"
#include "Gramm-Shmidt.hpp"

int
main()
{
    size_t m, n;
    //std::cin >> m >> n;
    //Matrix A(m, n);
    //std::cin >> A;
    std::vector<double> x = {1, -2, 3};
    std::vector<double> y = {1, 4, 10};
    std::vector<double> t = {4, 3, 0};
    std::vector<double> v = y;
    x = x / norm(x);
    y = y / norm(y);
    HouseholderMatrix H(y - x);
    std::cout << "prod = " << dot_product(H.w, v) << std::endl;
    matvec(H, v);
    for (size_t i = 0; i < v.size(); i++) {
        std::cout <<  v[i] / x[i] << ' ';
    }
    //for (auto i: z) std::cout << i << ' ';
    //std::cout << matrix_maximum_norm(A) << std::endl;
    //std::cout << norm(t) << std::endl;
    //std::cout << maximum_norm(x) << '\n';
    return 0;
}