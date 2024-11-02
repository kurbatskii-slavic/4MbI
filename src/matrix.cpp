#include "matrix.hpp"

double 
dot_product(std::vector<double>::const_iterator begin_x, std::vector<double>::const_iterator begin_y, size_t size) // vector dot_product
{
    double result = 0;
    for (size_t i = 0; i < size; i++) {
        result += (*begin_x) * (*begin_y);
        begin_x++;
        begin_y++;
    }
    return result;
}

double
norm(const std::vector<double> &x, size_t shift) // second Holder norm (||.||_2)
{
    return std::sqrt(dot_product(x.cbegin() + shift, x.cbegin() + shift, x.size() - shift));
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

std::ostream 
&operator<<(std::ostream &os, const std::vector<double>& v)
{
    for (auto i: v) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    return os;
}

std::ostream 
&operator<<(std::ostream &os, const HouseholderMatrix& T)
{
    std::cout << T.w;
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
operator+(const std::vector<double> &self, const std::vector<double> &other) // "x + y" overloading
{
    std::vector<double> res = self;
    res += other;
    return res;
}


void
operator+=(std::vector<double> &self, const std::vector<double> &other) // "x += y" overloading
{
    if (self.size() == other.size()) {
        for (size_t i = 0; i < other.size(); i++) {
            self[i] += other[i];
        }
    }
}

std::vector<double>
operator*(const Matrix &A, const std::vector<double> &x) // "A * v" overloading
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

Matrix 
operator*(const Matrix &A, const Matrix &B)
{
    Matrix C;
    for (size_t i = 0; i < B.cols; i++) {
        C[i] = A * B[i];
    }
    return C;
}

std::vector<double> 
operator*(const std::vector<double> x, const double &a) // "v * const" overloading
{
    std::vector<double> result = x;
    result *= a;
    return result;
}

void
operator*=(std::vector<double> &self, const double &a)
{
    for (size_t i = 0; i < self.size(); i++) {
        self[i] *= a;
    }   
}

std::vector<double> 
operator/(const std::vector<double> x, const double &a) // "v / const" overloading
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


std::vector<double> 
solve_triangular_system(const Matrix &R, const std::vector<double> &f) // Gauss method
{
    size_t n = f.size() - 1;
    std::vector<double> x(f.size());
    for (int m = n; m >= 0; m--) {
        double s = 0;
        for (size_t j = m + 1; j <= n; j++) {
            s += R(m, j) * x[j];
        }
        x[m] = (f[m] - s) / R(m, m);
    }
    return x;
}

void 
matvec(const HouseholderMatrix &H, std::vector<double> &v) // reflection operation (H * v)
{
    v -= H.w * 2 * dot_product(H.w.begin(), v.begin(), v.size());
    if (H.coeff == -1) {
        v *= H.coeff; 
    }
}

void 
matmul(const HouseholderMatrix &H, Matrix &A)
{
    for (size_t j = 0; j < A.cols; j++) {
        matvec(H, A[j]);
    }
}

void 
make_reflection(HouseholderMatrix &H, const std::vector<double> &x, const std::vector<double> &y, size_t shift) // Make H: x -> y
{
    double cos = dot_product(x.begin() + shift, y.begin() + shift, y.size()) / (norm(x, shift) * norm(y, shift));
    int sign = cos > 0 ? 1 : -1;
    std::vector<double> x_normalized = x;
    std::vector<double> y_normalized = y * sign;
    for (size_t i = 0; i < shift; i++) {
        x_normalized[i] = 0;
        y_normalized[i] = 0;
    }
    H.w = (x_normalized + y_normalized) / norm(x_normalized + y_normalized);
    H.coeff = -sign;
}

void
get_reflection(std::vector<double> &ref, std::vector<double> &x, size_t shift)
{
    double sub_vec_norm = norm(x, shift);
    for (size_t i = 0; i < shift; i++) {
        ref[i] = x[i];
    } 
    ref[shift] = sub_vec_norm;
}