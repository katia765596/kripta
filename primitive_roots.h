#ifndef PRIMITIVE_ROOTS_H
#define PRIMITIVE_ROOTS_H
#include <cstdint>
#include <vector>
class primitive_roots
{
public:
    static bool exists(uint64_t n);
    static std::vector<uint64_t> get_all(uint64_t n);
private:
    static uint64_t gcd(uint64_t a, uint64_t b);
    static uint64_t phi(uint64_t n);
    static uint64_t multiply_mod(uint64_t a, uint64_t b, uint64_t mod);
    static uint64_t power_mod(uint64_t a, uint64_t e, uint64_t mod);
    static std::vector<uint64_t> factorize(uint64_t n);//разложение n на простые (возвр вектор уникальных простых делителей)
};
#endif
