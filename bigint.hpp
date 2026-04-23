#ifndef BIGINT_HPP
#define BIGINT_HPP
#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>
class BigInt {
public:
    BigInt();
    BigInt(int64_t value);
    BigInt(const std::string& str);
    BigInt(const BigInt& other) = default;
    BigInt operator+(const BigInt& other) const;
    BigInt operator-(const BigInt& other) const;
    BigInt operator*(const BigInt& other) const;
    BigInt operator/(const BigInt& other) const;
    BigInt operator%(const BigInt& other) const;
    BigInt operator-() const;
    BigInt& operator+=(const BigInt& other);
    BigInt& operator-=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);
    BigInt& operator/=(const BigInt& other);
    BigInt& operator%=(const BigInt& other);
    bool operator==(const BigInt& other) const;
    bool operator!=(const BigInt& other) const;
    bool operator<(const BigInt& other) const;
    bool operator<=(const BigInt& other) const;
    bool operator>(const BigInt& other) const;
    bool operator>=(const BigInt& other) const;
    BigInt& operator++();
    BigInt operator++(int);
    BigInt& operator--();
    BigInt operator--(int);
    BigInt operator<<(size_t bits) const;
    BigInt operator>>(size_t bits) const;
    BigInt& operator<<=(size_t bits);
    BigInt& operator>>=(size_t bits);
    std::string to_string() const;
    int64_t to_int64() const;
    uint64_t to_uint64() const;
    bool is_zero() const;
    bool is_odd() const;
    bool is_negative() const;
    int sign() const;
    BigInt abs() const;
    static std::pair<BigInt, BigInt> divmod(const BigInt& a, const BigInt& b);
    static BigInt mod_pow(BigInt base, BigInt exp, const BigInt& mod);
    static BigInt mod_inverse(const BigInt& a, const BigInt& mod);
    static BigInt random(const BigInt& min, const BigInt& max);
    static bool is_prime(const BigInt& n, int certainty = 10);
    static BigInt generate_prime(int bits);
    static BigInt center_lift(const BigInt& x, const BigInt& mod);
private:
    std::vector<uint32_t> digits_;
    bool negative_;

    void trim();
    void from_string(const std::string& str);
    static int compare_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
    static std::vector<uint32_t> add_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
    static std::vector<uint32_t> sub_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
    static std::vector<uint32_t> mul_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
    static std::pair<std::vector<uint32_t>, std::vector<uint32_t>> divmod_abs(
        const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
};
#endif 