#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#ifndef MATRIX_HPP
#define MATRIX_HPP

struct HouseholderMatrix // struct for Householder matrices
{
    std::vector<double> w;
    size_t n = 0;
    HouseholderMatrix(const std::vector<double> &v): w(v / norm(v)), n(w.size()) {} // constructor

};


struct Matrix // struct for matrices
{
    std::vector<double> arr; // array of columns (convenient for calculations)
    const size_t rows, cols; // dimensions

    Matrix(size_t m, size_t n) : rows(m), cols(n), arr(std::vector<double>(m * n, 0)) {} // constructors
    Matrix(size_t n) : Matrix(n, n) {}  // square matrix constructor
    double &operator()(size_t i, size_t j) { return arr[i * cols + j]; } // element access
    Matrix operator-(Matrix B) const {
        Matrix C = *this;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                C(i, j) -= B(i, j);
            }    
        }
        return C;
    }
};

double dot_product(const std::vector<double> &x, const std::vector<double> &y);
double norm(const std::vector<double> &x);
std::vector<double> operator*(Matrix &A, const std::vector<double> &x);
std::vector<double> operator-(const std::vector<double> &self, const std::vector<double> &other);
std::vector<double> operator*(const std::vector<double> x, double a);
std::vector<double> operator/(const std::vector<double> x, double a);
std::ostream &operator<<(std::ostream &os, Matrix& A);
std::istream &operator>>(std::istream &is, Matrix& A);
double maximum_norm(std::vector<double> x);

#endif
