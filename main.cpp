#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>

#include "matrix.hpp"

int
main()
{
    Matrix A(3);
    std::cin >> A;
    std::vector<double> x = {1, -2, 3};
    std::vector<double> y = {1, -4, 10}, z = A * x;
    for (auto i: z) std::cout << i << ' ';
    std::cout << std::endl;
    std::cout << maximum_norm(x) << '\n';
    return 0;
}