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
    std::vector<double> x = {1, -2, 3, 4};
    std::cout << maximum_norm(x) << '\n';
    return 0;
}