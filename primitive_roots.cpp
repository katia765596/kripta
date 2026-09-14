#include "primitive_roots.h"
#include <stdexcept>
uint64_t primitive_roots::gcd(uint64_t a, uint64_t b)
{
    while (b != 0)
    {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}
uint64_t primitive_roots::phi(uint64_t n)
{
    uint64_t result = n;
    for (uint64_t p : factorize(n))//перебирает все простые множители числа n
    {
        result = result / p * (p - 1);//ф эйлера ф(n)=n*П(1-1/p) для каждого простого р
    }
    return result;
}
uint64_t primitive_roots::multiply_mod(uint64_t a, uint64_t b, uint64_t mod)
{
    uint64_t result = 0;
    a %= mod;
    while (b != 0)
    {
        if (b & 1)
            result = result >= mod - a ? result - (mod - a) : result + a;
        a = a >= mod - a ? a - (mod - a) : a + a;
        b >>= 1;
    }
    return result;
}
uint64_t primitive_roots::power_mod(uint64_t a, uint64_t e, uint64_t mod)
{
    uint64_t result = 1 % mod;
    a %= mod;
    while (e != 0)
    {
        if (e & 1)
            result = multiply_mod(result, a, mod);
        a = multiply_mod(a, a, mod);
        e >>= 1;
    }
    return result;
}
std::vector<uint64_t> primitive_roots::factorize(uint64_t n)//разлагает n на простые 
{
    std::vector<uint64_t> result;
    for (uint64_t p = 2; p <= n / p; ++p)//усл экв р*р<=n
    {
        if (n % p == 0)
        {
            result.push_back(p);
            while (n % p == 0)
                n /= p;
        }
    }
    if (n > 1)
        result.push_back(n);
    return result;
}
bool primitive_roots::exists(uint64_t n)
{
    if (n == 2 || n == 4)
        return true;
    if (n < 2)
        return false;
    uint64_t odd_part = n;
    size_t twos = 0;
    while ((odd_part & 1) == 0)//пока n четное
    {
        odd_part >>= 1;
        ++twos;
    }
    if (twos > 1)
        return false;
    if (odd_part == 1)//если после выделения всех двоек осталось 1, n=2^k
        return false;
    return factorize(odd_part).size() == 1;
}//для существования первообразного корня необходимо, чтобы нечетная часть была степенью нечетного простого числа (т е имела ровно один простой множитель)
std::vector<uint64_t> primitive_roots::get_all(uint64_t n)
{
    if (!exists(n))
        return {};
    if (n == 2)
        return { 1 };
    if (n == 4)
        return { 3 };
    uint64_t order = phi(n);//порядок мультипликативной группы по модулю n равен фи
    std::vector<uint64_t> factors = factorize(order);//разлагаем порядок на простые
    uint64_t generator = 0;
    for (uint64_t g = 2; g < n; ++g)
    {
        if (gcd(g, n) != 1)
            continue;
        bool good = true;
        for (uint64_t q : factors)//для каждого простого делителя q порчдка, проверяем что g^(order/q!=1(mod n)
        {
            if (power_mod(g, order / q, n) == 1)
            {
                good = false;
                break;
            }
        }
        if (good)
        {
            generator = g;
            break;
        }
    }
    if (generator == 0)
        throw std::runtime_error("primitive root not found");
    std::vector<uint64_t> result;//первообразные имеют вид generator^i,где i взаимно просто с order
    for (uint64_t i = 1; i <= order; ++i)
    {
        if (gcd(i, order) == 1)
            result.push_back(power_mod(generator, i, n));
    }
    return result;
}
