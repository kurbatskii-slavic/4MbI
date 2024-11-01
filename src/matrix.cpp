#include "matrix.hpp"


double 
dot_product(const std::vector<double> &x, const std::vector<double> &y)
{
    double result = 0;
    if (x.size() == y.size()) {
        for (size_t i = 0; i < y.size(); i++) {
            result += x[i] * y[i];
        }
    }
    return result;
}

double
norm(const std::vector<double> &x)
{
    return std::sqrt(dot_product(x, x));
}


std::ostream 
&operator<<(std::ostream &os, Matrix& A) // cout overloading
{
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            os << A(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

std::istream 
&operator>>(std::istream &is, Matrix& A) // cin overloading
{
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            is >> A(i, j);
        }
    }
    return is;
}


std::vector<double> 
operator-(const std::vector<double> &self, const std::vector<double> &other)
{
    std::vector<double> res(other.size(), 0);
    if (self.size() == other.size()) {
        for (size_t i = 0; i < other.size(); i++) {
            res[i] = self[i] - other[i];
        }
    }
    return res;
}


std::vector<double>
operator*(Matrix &A, const std::vector<double> &x)
{
    std::vector<double> result(x.size(), 0);
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0;
        for (size_t j = 0; j < A.cols; j++) {
            sum += x[j] * A(i, j);           
        }    
        result[i] = sum;
    }
    return result;
}

std::vector<double> 
operator*(const std::vector<double> x, double a)
{
    std::vector<double> result = x;
    for (size_t i = 0; i < x.size(); i++) {
        result[i] *= a;
    }
    return result;
}
std::vector<double> 
operator/(const std::vector<double> x, double a)
{
    std::vector<double> result = x;
    for (size_t i = 0; i < x.size(); i++) {
        result[i] /= a;
    }
    return result;
}


double
maximum_norm(std::vector<double> x)
{
    double m = 0;    
    for (auto x_i: x) {
        double abs = x_i > 0 ? x_i : -x_i;
        m = abs > m ? abs : m;
    }
    return m;
}