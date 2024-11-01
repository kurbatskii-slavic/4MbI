#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <chrono>

#include "matrix.hpp"

int
main()
{
    Matrix A(3);
    //std::cin >> A;
    std::vector<double> x = {1, -2, 3};
    std::vector<double> y = {1, 4, 10}, z = x;
    std::vector<double> t = {4, 3, 0};
    for (auto i: z) std::cout << i << ' ';
    //std::cout << dot_product(x, y) << std::endl;
    //std::cout << norm(t) << std::endl;
    //std::cout << maximum_norm(x) << '\n';
    return 0;
}