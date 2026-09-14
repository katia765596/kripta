#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include "big_int.h"
#include <string>
#include <vector>
class polynomial
{
public:
    using coeff_type = big_int;
    using container = std::vector<coeff_type>;
    polynomial();
    explicit polynomial(
        size_t degree,
        const coeff_type& value = 0);
    explicit polynomial(const container& coeffs);
    size_t degree() const;
    const coeff_type& operator[](size_t i) const;//для чтения
    coeff_type& operator[](size_t i);//возвр ссылку на коэф по индексу i, для записи
    polynomial operator+(const polynomial& other) const;
    polynomial operator-(const polynomial& other) const;
    polynomial operator*(const polynomial& other) const;
    polynomial operator*(const coeff_type& scalar) const;
    polynomial mod(const coeff_type& mod) const;
    polynomial center_lift(const coeff_type& mod) const;
    polynomial mul_mod(
        const polynomial& other,
        const coeff_type& mod,
        size_t N) const;
    bool operator==(const polynomial& other) const;
    bool operator!=(const polynomial& other) const;
    polynomial inverse_mod(
        const coeff_type& mod,
        size_t N) const;
    std::string to_string() const;
    static polynomial random(
        size_t degree,
        const std::vector<coeff_type>& coeff_set);
    void resize(
        size_t N,
        const coeff_type& fill = 0);
    container get_coeffs() const
    {
        return coeffs_;
    }
    void set_coeffs(const container& coeffs)
    {
        coeffs_ = coeffs;
        trim();
    }
private:
    container coeffs_;
    void trim();
};
#endif
