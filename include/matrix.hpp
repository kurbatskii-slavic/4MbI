#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#ifndef MATRIX_HPP
#define MATRIX_HPP


struct Matrix;
struct HouseholderMatrix;
double norm(const std::vector<double> &x, size_t shift=0);
std::vector<double> operator*(Matrix &A, const std::vector<double> &x);
std::vector<double> operator-(const std::vector<double> &self, const std::vector<double> &other);
void operator-=(std::vector<double> &self, const std::vector<double> &other);
std::vector<double> operator+(const std::vector<double> &self, const std::vector<double> &other);
void operator+=(std::vector<double> &self, const std::vector<double> &other);
std::vector<double> operator*(const std::vector<double> x, const double &a);
void operator*=(std::vector<double> &self, const double &a);
std::vector<double> operator/(const std::vector<double> x, const double &a);

struct HouseholderMatrix // struct for Householder matrices
{
    std::vector<double> w; // normal vector
    int coeff = 1;
    HouseholderMatrix(){};
    HouseholderMatrix(const std::vector<double> &v): w(v / norm(v)) {} // constructor
};


struct Matrix // struct for matrices
{
    std::vector<std::vector<double>> arr; // array of columns (convenient for calculations)
    size_t rows, cols; // dimensions

    Matrix(size_t m, size_t n) : rows(m), cols(n), arr(n, std::vector<double>(m, 0)) {} // constructors
    double &operator()(size_t i, size_t j) { return arr[j][i]; } // element access
    double operator()(size_t i, size_t j) const { return arr[j][i]; } // element access
    Matrix operator-(Matrix B) const {
        Matrix C = *this;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                C(i, j) -= B(i, j);
            }    
        }
        return C;
    }
    std::vector<double> &operator[](size_t i) { return arr[i]; } // get i-column
};


double dot_product(std::vector<double>::const_iterator begin_x, std::vector<double>::const_iterator begin_y, size_t size);
std::ostream &operator<<(std::ostream &os, const Matrix& A);
std::ostream &operator<<(std::ostream &os, const std::vector<double>& v);
std::ostream &operator<<(std::ostream &os, const HouseholderMatrix& T);
std::istream &operator>>(std::istream &is, Matrix& A);
double maximum_norm(const std::vector<double> &x);
double matrix_maximum_norm(const Matrix &A);

void matvec(const HouseholderMatrix &H, std::vector<double> &v);
void make_reflection(HouseholderMatrix &H, const std::vector<double> &x, const std::vector<double> &y, size_t shift);
void get_reflection(std::vector<double> &ref, std::vector<double> &x, size_t shift);

#endif
