#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cassert>
#include <iostream>

namespace primitive_roots_ns {
    using uint64 = std::uint64_t;
    uint64 gcd(uint64 a, uint64 b) {
        while (b != 0) {
            uint64 t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
    uint64 power_mod(uint64 base, uint64 exp, uint64 mod) {
        if (mod == 1) return 0;
        uint64 result = 1;
        base %= mod;
        while (exp > 0) {
            if (exp & 1) {
                result = (uint64)((__int128)result * base % mod);
            }
            base = (uint64)((__int128)base * base % mod);
            exp >>= 1;
        }
        return result;
    }
    bool is_prime(uint64 n) {
        if (n < 2) return false;
        if (n == 2 || n == 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        uint64 limit = static_cast<uint64>(std::sqrt(n)) + 1;
        for (uint64 i = 5; i <= limit; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0) return false;
        }
        return true;
    }
    std::vector<uint64> factorize(uint64 n) {
        std::vector<uint64> factors;
        if (n == 1) return factors;
        while (n % 2 == 0) {
            factors.push_back(2);
            n /= 2;
        }
        uint64 d = 3;
        uint64 limit = static_cast<uint64>(std::sqrt(n)) + 1;
        while (d <= limit && n > 1) {
            while (n % d == 0) {
                factors.push_back(d);
                n /= d;
            }
            d += 2;
            limit = static_cast<uint64>(std::sqrt(n)) + 1;
        }
        if (n > 1) factors.push_back(n);
        return factors;
    }
    std::vector<uint64> unique_prime_factors(uint64 n) {
        auto factors = factorize(n);
        std::sort(factors.begin(), factors.end());
        factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
        return factors;
    }
    uint64 euler_phi(uint64 n) {
        uint64 result = n;
        auto primes = unique_prime_factors(n);
        for (uint64 p : primes) {
            result = result / p * (p - 1);
        }
        return result;
    }
    bool exists_primitive_root(uint64 n) {
        if (n == 2 || n == 4) return true;
        auto factors = factorize(n);
        if (factors.empty()) return false;
        std::sort(factors.begin(), factors.end());
        uint64 distinct = 1;
        for (size_t i = 1; i < factors.size(); ++i) {
            if (factors[i] != factors[i - 1]) ++distinct;
        }
        if (distinct == 1) {
            uint64 p = factors[0];
            if (p == 2) return false; 
            return true;
        }
        if (distinct == 2) {
            if (factors[0] != 2) return false;
            if (factors.size() > 1 && factors[1] == 2) return false;
            return true;
        }
        return false;
    }
    uint64 find_one_primitive_root(uint64 n) {
        uint64 phi = euler_phi(n);
        auto prime_factors_of_phi = unique_prime_factors(phi);
        for (uint64 a = 2; a < n; ++a) {
            if (gcd(a, n) != 1) continue;
            bool ok = true;
            for (uint64 q : prime_factors_of_phi) {
                if (power_mod(a, phi / q, n) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return a;
        }
        return 0; 
    }
    std::vector<uint64> primitive_roots(uint64 n) {
        if (n < 2) return {};
        if (!exists_primitive_root(n)) return {};
        if (n == 2) return { 1 };
        uint64 phi = euler_phi(n);
        uint64 g = find_one_primitive_root(n);
        if (g == 0) return {};
        std::vector<uint64> roots;
        uint64 current = 1;
        for (uint64 i = 1; i <= phi; ++i) {
            current = (uint64)((__int128)current * g % n);
            if (gcd(i, phi) == 1) {
                roots.push_back(current);
            }
        }
        std::sort(roots.begin(), roots.end());
        return roots;
    }

}
void test_primitive_roots() {
    using namespace primitive_roots_ns;
    auto test = [](uint64 n, const std::vector<uint64>& expected) {
        auto got = primitive_roots(n);
        if (got == expected) {
            std::cout << "PASS: n = " << n << std::endl;
        }
        else {
            std::cout << "FAIL: n = " << n << "\n  Expected: ";
            for (auto x : expected) std::cout << x << " ";
            std::cout << "\n  Got: ";
            for (auto x : got) std::cout << x << " ";
            std::cout << std::endl;
            std::abort();
        }
        };
    test(2, { 1 });
    test(3, { 2 });
    test(4, { 3 });
    test(5, { 2, 3 });
    test(6, { 5 });
    test(7, { 3, 5 });
    test(9, { 2, 5 });
    test(10, { 3, 7 });
    test(14, { 3, 5 });
    test(18, { 5, 11 });
    test(25, { 2, 3, 8, 12, 13, 17, 22, 23 });     
    test(27, { 2, 5, 11, 14, 20, 23 });          
    test(49, { 3, 5, 10, 12, 17, 24, 26, 33, 38, 40, 45, 47 });  
    test(50, { 3, 13, 17, 23, 27, 33, 37, 47 });                 
    test(1, {});
    test(8, {});
    test(12, {});
    test(15, {});
    test(16, {});
    test(20, {});
    test(21, {});
    test(24, {});
    test(28, {});
    test(30, {});
    test(36, {});
    std::cout << "\nAll tests passed successfully!\n";
}
int main() {
    test_primitive_roots();
    return 0;
}