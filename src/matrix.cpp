#include "matrix.hpp"

double 
dot_product(const std::vector<double> &x, const std::vector<double> &y) // vector dot_product
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
norm(const std::vector<double> &x) // second Holder norm (||.||_2)
{
    return std::sqrt(dot_product(x, x));
}

std::ostream 
&operator<<(std::ostream &os, const Matrix& A) // cout overloading
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
operator-(const std::vector<double> &self, const std::vector<double> &other) // "x - y" overloading
{
    std::vector<double> res = self;
    res -= other;
    return res;
}


void
operator-=(std::vector<double> &self, const std::vector<double> &other) // "x -= y" overloading
{
    if (self.size() == other.size()) {
        for (size_t i = 0; i < other.size(); i++) {
            self[i] -= other[i];
        }
    }
}

std::vector<double>
operator*(Matrix &A, const std::vector<double> &x) // "A * v" overloading
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
operator*(const std::vector<double> x, double a) // "v * const" overloading
{
    std::vector<double> result = x;
    for (size_t i = 0; i < x.size(); i++) {
        result[i] *= a;
    }
    return result;
}

std::vector<double> 
operator/(const std::vector<double> x, double a) // "v / const" overloading
{
    std::vector<double> result = x;
    for (size_t i = 0; i < x.size(); i++) {
        result[i] /= a;
    }
    return result;
}

double
maximum_norm(const std::vector<double> &x) // vector maximum_norm
{
    double m = 0;    
    for (auto x_i: x) {
        double abs = x_i > 0 ? x_i : -x_i;
        m = abs > m ? abs : m;
    }
    return m;
}

double
matrix_maximum_norm(const Matrix &A) // matrix maximum norm
{
    double m = 0;
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0;
        for (size_t j = 0; j < A.cols; j++) {
            double abs = A(i, j) > 0 ? A(i, j) : -A(i, j);
            sum += abs;
        }
        m = sum > m ? sum : m;
    } 
    return m;
}

void 
matvec(const HouseholderMatrix &H, std::vector<double> &v) // reflection operation (H * v)
{
    v -= H.w * 2 * dot_product(H.w, v); 
}