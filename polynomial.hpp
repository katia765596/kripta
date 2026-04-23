#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP
#include "bigint.hpp"
#include <vector>
#include <string>
class Polynomial {
public:
    using CoeffType = BigInt;
    using Container = std::vector<CoeffType>;
    Polynomial();
    explicit Polynomial(size_t degree, const CoeffType& value = 0);
    explicit Polynomial(const Container& coeffs);
    size_t degree() const;
    const CoeffType& operator[](size_t i) const;
    CoeffType& operator[](size_t i);
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;
    Polynomial operator*(const CoeffType& scalar) const;
    Polynomial mod(const CoeffType& mod) const;
    Polynomial center_lift(const CoeffType& mod) const;
    Polynomial mul_mod(const Polynomial& other, const CoeffType& mod, size_t N) const;
    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;
    Polynomial inverse_mod(const CoeffType& mod, size_t N) const;
    std::string to_string() const;
    static Polynomial random(size_t degree, const std::vector<CoeffType>& coeff_set);
    void resize(size_t N, const CoeffType& fill = 0);
    Container get_coeffs() const { return coeffs_; }
    void set_coeffs(const Container& new_coeffs) { coeffs_ = new_coeffs; trim(); }
private:
    Container coeffs_;
    void trim();
};
#endif 