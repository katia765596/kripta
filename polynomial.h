#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
class polynomial//полиномы над полем gf(2) набор 64 битных слов, где каждый бит соотв коэф при соотв степени
{
public:
    polynomial();
    explicit polynomial(uint64_t value);
    explicit polynomial(const std::vector<uint64_t>& value);//констр из вектора 64битных слов
    static polynomial from_binary(const std::string& value);//полином из строки двоичного вида
    static polynomial one();
    static polynomial x();
    bool is_zero() const;
    bool is_one() const;
    size_t degree() const;
    bool get_bit(size_t index) const;//возвр знач бита
    void set_bit(size_t index, bool value);//устанавл бит на поз индекс в знач value
    uint64_t to_uint64() const;//преобраз полином в 64битное целое
    std::string to_binary() const;
    polynomial add(const polynomial& other) const;
    polynomial multiply(const polynomial& other) const;
    polynomial mod(const polynomial& modulus) const;
    polynomial multiply_mod(const polynomial& other, const polynomial& modulus) const;
    polynomial inverse(const polynomial& modulus) const;
    polynomial power_mod(uint64_t exponent, const polynomial& modulus) const;
    polynomial gcd(const polynomial& other) const;
    bool is_irreducible() const;//проверка на неприводимость
    bool operator==(const polynomial& other) const;
    bool operator!=(const polynomial& other) const;
private:
    std::vector<uint64_t> data;
    void normalize();//удаляет старшие нулевые биты
    void xor_shifted(const polynomial& other, size_t shift);
    polynomial divide_remainder(const polynomial& divisor, polynomial& remainder) const;
};
class finite_field//представляет конечное поле gf(2^m)
{
public:
    explicit finite_field(const polynomial& modulus);//неприводимый полином степени m сохраняет его как модуль
    polynomial add(const polynomial& first, const polynomial& second) const;
    polynomial multiply(const polynomial& first, const polynomial& second) const;
    polynomial inverse(const polynomial& value) const;
    polynomial power(const polynomial& value, uint64_t exponent) const;
    const polynomial& get_modulus() const;//возвр конст ссылку на модуль поля
private:
    polynomial modulus;
};
#endif
