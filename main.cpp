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
    size_t n;
    std::cin >> n;
    Matrix A(n, n);
    std::cin >> A;
    std::cout << '\n';
    std::vector<HouseholderMatrix> H;
   // std::cout << std::endl << A << std::endl;
    QRHouseholder(H, A);
    std::pair<std::vector<HouseholderMatrix>, Matrix> result = std::make_pair(H, A);
    std::cout << result << std::endl;
    std::reverse(H.begin(), H.end());
    for (auto T: H) {
        for (size_t j = 0; j < A.cols; j++) {
            matvec(T, A[j]);
        }
    }
    std::cout << A;
//    vector<double> x = {-4 ,1, 4, 3, 1};
  //  std::cout << x << std::endl;
    
    return 0;
}