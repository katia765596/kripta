#include "big_int.h"
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/random.hpp>
#include <chrono>
#include <stdexcept>
namespace
{
    boost::random::mt19937& rng()//функция возвращает ссылку на статический объект генератора
    {
        static boost::random::mt19937 generator(
            static_cast<unsigned>(std::chrono::steady_clock::now().time_since_epoch().count()));//объект иниц 1раз при первом вызове, аргемуент - тек время в наносекундах
        return generator;
    }
}
big_int random_bigint(const big_int& min, const big_int& max)
{
    if (min > max) {throw std::invalid_argument("random_bigint: min > max");}
    boost::random::uniform_int_distribution<big_int> distribution(min, max);
    return distribution(rng());
}
big_int mod_pow(const big_int& base, const big_int& exp, const big_int& mod)
{
    if (mod <= 0) {throw std::invalid_argument("mod_pow: mod must be positive");}
    return boost::multiprecision::powm(base, exp, mod);
}
big_int mod_inverse(const big_int& a, const big_int& mod)
{
    if (mod <= 1) {throw std::invalid_argument("mod_inverse: invalid modulus");}
    big_int t = 0;
    big_int new_t = 1;
    big_int r = mod;
    big_int new_r = a % mod;
    if (new_r < 0)
    {new_r += mod;}
    while (new_r != 0)
    {
        big_int q = r / new_r;
        big_int tmp = t;
        t = new_t;
        new_t = tmp - q * new_t;
        tmp = r;
        r = new_r;
        new_r = tmp - q * new_r;
    }
    if (r != 1) {throw std::runtime_error("mod_inverse: not invertible");}
    t %= mod;
    if (t < 0)
    {t += mod;}
    return t;
}
bool is_prime(const big_int& n, int certainty)
{
    if (n < 2)
    {return false;}
    if (n == 2 || n == 3)
    {return true;}
    if (n % 2 == 0)
    {return false;}
    const int small_primes[] = {
        3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71
    };
    for (int p : small_primes)
    {
        if (n == p)
        {return true;}
        if (n % p == 0)
        {return false;}
    }
    return boost::multiprecision::miller_rabin_test(n, certainty);
}
big_int generate_prime(int bits)
{
    if (bits < 2)
    {throw std::invalid_argument("bits must be >= 2");}
    const big_int min = big_int(1) << (bits - 1);
    const big_int max = (big_int(1) << bits) - 1;
    while (true)
    {
        big_int candidate = random_bigint(min, max) | 1;
        if (is_prime(candidate, 25))
        {return candidate;}
    }
}
big_int center_lift(const big_int& x, const big_int& mod)
{
    big_int r = x % mod;
    if (r < 0)
    {r += mod;}
    const big_int half = mod / 2;
    if (r > half)
    {r -= mod;}
    return r;//центрированное значение x в диапозоне -mod/2 до mod/2
}
