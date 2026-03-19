#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <tuple>
#include <cmath>
#include <stdexcept>
#include <algorithm>
class MathService {
public:
    static int legendreSymbol(long long a, long long p) {
        if (p <= 2 || a < 0) throw std::invalid_argument("Invalid input");
        a %= p;
        if (a == 0) return 0;
        if (a == 1) return 1;

        if (a % 2 == 0) {
            int s = legendreSymbol(a / 2, p);
            int r = (p * p - 1) / 8;
            return (r % 2 == 0 ? s : -s);
        }
        int res = legendreSymbol(p % a, a);
        if (((a - 1) * (p - 1) / 4) % 2) res = -res;
        return res;
    }
    static int jacobiSymbol(long long a, long long n) {
        if (n <= 0 || n % 2 == 0) throw std::invalid_argument("Invalid input");
        a %= n;
        int result = 1;
        while (a != 0) {
            while (a % 2 == 0) {
                a /= 2;
                if ((n % 8 == 3) || (n % 8 == 5)) result = -result;
            }
            std::swap(a, n);
            if ((a % 4 == 3) && (n % 4 == 3)) result = -result;
            a %= n;
        }
        return (n == 1 ? result : 0);
    }
    static long long gcd(long long a, long long b) {
        while (b != 0) {
            long long t = b;
            b = a % b;
            a = t;
        }
        return std::abs(a);
    }
    static std::tuple<long long, long long, long long> extendedGCD(long long a, long long b) {
        long long x0 = 1, y0 = 0;
        long long x1 = 0, y1 = 1;
        while (b != 0) {
            long long q = a / b;
            long long r = a % b;
            a = b;
            b = r;
            long long x = x0 - q * x1;
            long long y = y0 - q * y1;

            x0 = x1; x1 = x;
            y0 = y1; y1 = y;
        }
        return { a, x0, y0 };
    }
    static long long modPow(long long base, long long exp, long long mod) {
        if (mod <= 0) throw std::invalid_argument("Invalid modulus");
        base %= mod;
        long long result = 1;
        while (exp > 0) {
            if (exp & 1) result = (result * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }
        return result;
    }
    static long long modInverse(long long a, long long mod) {
        auto [g, x, y] = extendedGCD(a, mod);
        if (g != 1) throw std::invalid_argument("Inverse does not exist");
        return (x % mod + mod) % mod;
    }
    static long long eulerPhi(long long n) {
        long long count = 0;
        for (long long i = 1; i <= n; ++i) {
            if (gcd(i, n) == 1) ++count;
        }
        return count;
    }
    static long long eulerPhiFactor(long long n) {
        long long result = n;
        for (long long p = 2; p * p <= n; ++p) {
            if (n % p == 0) {
                while (n % p == 0) n /= p;
                result -= result / p;
            }
        }
        if (n > 1) result -= result / n;
        return result;
    }
    static double eulerPhiDFT(long long n) {
        double sum = 0;
        for (long long k = 1; k <= n; ++k) {
            sum += gcd(k, n) * std::cos(2 * M_PI * k / n);
        }
        return sum;
    }
};
int main() {
    std::cout << "=== MathService Demo ===\n";
    std::cout << "GCD(48,180) = " << MathService::gcd(48, 180) << "\n";
    auto [g, x, y] = MathService::extendedGCD(48, 180);
    std::cout << "ExtendedGCD: gcd=" << g << " x=" << x << " y=" << y << "\n";
    std::cout << "ModPow(3,200,13) = " << MathService::modPow(3, 200, 13) << "\n";
    std::cout << "ModInverse(3,11) = " << MathService::modInverse(3, 11) << "\n";
    std::cout << "EulerPhi(36) = " << MathService::eulerPhi(36) << "\n";
    std::cout << "EulerPhiFactor(36) = " << MathService::eulerPhiFactor(36) << "\n";
    std::cout << "EulerPhiDFT(36) = " << MathService::eulerPhiDFT(36) << "\n";
    std::cout << "Legendre(5,11) = " << MathService::legendreSymbol(5, 11) << "\n";
    std::cout << "Jacobi(1001,9907) = " << MathService::jacobiSymbol(1001, 9907) << "\n";
    return 0;
}

