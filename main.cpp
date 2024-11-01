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
using std::vector;

int
main()
{
    size_t m, n;
    //std::cin >> m >> n;
    Matrix A(m, n);
    //std::cin >> A;
    HouseholderMatrix H;
    vector<double> x = {-1 , 1};
    vector<double> y = {1, 0};
    make_reflection(H, x, y);
    for(auto i: H.w) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    return 0;
}