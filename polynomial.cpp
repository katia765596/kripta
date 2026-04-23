#include "polynomial.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>
Polynomial::Polynomial() = default;
Polynomial::Polynomial(size_t degree, const CoeffType& value) : coeffs_(degree + 1, value) { trim(); }
Polynomial::Polynomial(const Container& coeffs) : coeffs_(coeffs) { trim(); }
void Polynomial::trim() {
    while (!coeffs_.empty() && coeffs_.back() == 0) coeffs_.pop_back();
}
size_t Polynomial::degree() const { return coeffs_.empty() ? 0 : coeffs_.size() - 1; }
const BigInt& Polynomial::operator[](size_t i) const { return coeffs_[i]; }
BigInt& Polynomial::operator[](size_t i) { return coeffs_[i]; }
Polynomial Polynomial::operator+(const Polynomial& other) const {
    size_t max_size = std::max(coeffs_.size(), other.coeffs_.size());
    Container res(max_size, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i) res[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.size(); ++i) res[i] += other.coeffs_[i];
    return Polynomial(res);
}
Polynomial Polynomial::operator-(const Polynomial& other) const {
    size_t max_size = std::max(coeffs_.size(), other.coeffs_.size());
    Container res(max_size, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i) res[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.size(); ++i) res[i] -= other.coeffs_[i];
    return Polynomial(res);
}
Polynomial Polynomial::operator*(const Polynomial& other) const {
    if (coeffs_.empty() || other.coeffs_.empty()) return Polynomial();
    Container res(coeffs_.size() + other.coeffs_.size() - 1, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i)
        for (size_t j = 0; j < other.coeffs_.size(); ++j)
            res[i + j] += coeffs_[i] * other.coeffs_[j];
    return Polynomial(res);
}
Polynomial Polynomial::operator*(const CoeffType& scalar) const {
    if (scalar == 0) return Polynomial();
    Container res(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i) res[i] = coeffs_[i] * scalar;
    return Polynomial(res);
}
Polynomial Polynomial::mod(const CoeffType& mod) const {
    Container res(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i) {
        res[i] = coeffs_[i] % mod;
        if (res[i] < 0) res[i] += mod;
    }
    return Polynomial(res);
}
Polynomial Polynomial::center_lift(const CoeffType& mod) const {
    Container res(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
        res[i] = BigInt::center_lift(coeffs_[i], mod);
    return Polynomial(res);
}
Polynomial Polynomial::mul_mod(const Polynomial& other, const CoeffType& mod, size_t N) const {
    Container res(N, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i) {
        if (coeffs_[i] == 0) continue;
        for (size_t j = 0; j < other.coeffs_.size(); ++j) {
            if (other.coeffs_[j] == 0) continue;
            size_t idx = (i + j) % N;
            BigInt prod = coeffs_[i] * other.coeffs_[j];
            BigInt add = (res[idx] + prod) % mod;
            if (add < 0) add += mod;
            res[idx] = add;
        }
    }
    return Polynomial(res);
}
bool Polynomial::operator==(const Polynomial& other) const { return coeffs_ == other.coeffs_; }
bool Polynomial::operator!=(const Polynomial& other) const { return !(*this == other); }
std::string Polynomial::to_string() const {
    std::ostringstream oss;
    for (size_t i = 0; i < coeffs_.size(); ++i) {
        if (i > 0) oss << " + ";
        oss << coeffs_[i].to_string() << "x^" << i;
    }
    if (coeffs_.empty()) oss << "0";
    return oss.str();
}
Polynomial Polynomial::random(size_t degree, const std::vector<CoeffType>& coeff_set) {
    if (coeff_set.empty()) throw std::invalid_argument("empty coeff set");
    Container coeffs(degree + 1);
    for (size_t i = 0; i <= degree; ++i) {
        size_t idx = BigInt::random(0, coeff_set.size() - 1).to_uint64();
        coeffs[i] = coeff_set[idx];
    }
    return Polynomial(coeffs);
}
void Polynomial::resize(size_t N, const CoeffType& fill) {
    coeffs_.resize(N, fill);
    trim();
}
Polynomial Polynomial::inverse_mod(const CoeffType& mod, size_t N) const {
    Polynomial a = this->mod(mod);
    a.resize(N, 0);
    std::vector<std::vector<BigInt>> mat(N, std::vector<BigInt>(N));
    for (size_t k = 0; k < N; ++k)
        for (size_t i = 0; i < N; ++i)
            mat[k][i] = a.coeffs_[(k + N - i) % N];
    std::vector<BigInt> b(N, 0);
    b[0] = 1;
    for (size_t col = 0; col < N; ++col) {
        size_t pivot = col;
        while (pivot < N && mat[pivot][col] == 0) ++pivot;
        if (pivot == N) throw std::runtime_error("Matrix singular");
        if (pivot != col) {
            std::swap(mat[col], mat[pivot]);
            std::swap(b[col], b[pivot]);
        }
        BigInt inv_pivot = BigInt::mod_inverse(mat[col][col], mod);
        for (size_t j = col; j < N; ++j) {
            mat[col][j] = (mat[col][j] * inv_pivot) % mod;
            if (mat[col][j] < 0) mat[col][j] += mod;
        }
        b[col] = (b[col] * inv_pivot) % mod;
        if (b[col] < 0) b[col] += mod;
        for (size_t row = 0; row < N; ++row) {
            if (row != col && mat[row][col] != 0) {
                BigInt factor = mat[row][col];
                for (size_t j = col; j < N; ++j) {
                    mat[row][j] = (mat[row][j] - factor * mat[col][j]) % mod;
                    if (mat[row][j] < 0) mat[row][j] += mod;
                }
                b[row] = (b[row] - factor * b[col]) % mod;
                if (b[row] < 0) b[row] += mod;
            }
        }
    }
    Container inv_coeffs(N);
    for (size_t i = 0; i < N; ++i) {
        inv_coeffs[i] = b[i] % mod;
        if (inv_coeffs[i] < 0) inv_coeffs[i] += mod;
    }
    Polynomial inv(inv_coeffs);
    inv.trim();
    return inv;
}