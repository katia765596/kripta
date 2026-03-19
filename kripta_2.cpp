#include <iostream>
#include <tuple>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <cassert>

class MathService {
public:
    static long long modPow(long long base, long long exponent, long long mod) {
        if (mod <= 0) throw std::invalid_argument("Modulo must be positive");
        base %= mod;
        long long result = 1;
        while (exponent > 0) {
            if (exponent & 1) result = (result * base) % mod;
            base = (base * base) % mod;
            exponent >>= 1;
        }
        return result;
    }
};

class IProbabilisticPrimalityTest {
public:
    virtual ~IProbabilisticPrimalityTest() = default;
    virtual bool isPrime(long long n, double minProbability) = 0;
};
class ProbabilisticPrimalityTest : public IProbabilisticPrimalityTest {
protected:
    int computeIterations(double minProbability) const {
        if (minProbability < 0.5 || minProbability >= 1.0)
            throw std::invalid_argument("Probability must be in [0.5, 1)");
        return static_cast<int>(std::ceil(std::log(1.0 - minProbability) / std::log(0.5)));
    }
    virtual bool singleIteration(long long n) = 0;
public:
    bool isPrime(long long n, double minProbability) override {
        if (n < 2) return false;
        int iterations = computeIterations(minProbability);
        for (int i = 0; i < iterations; ++i) {
            if (!singleIteration(n)) return false;
        }
        return true;
    }
};
class FermatPrimalityTest : public ProbabilisticPrimalityTest {
private:
    std::mt19937_64 rng;

protected:
    bool singleIteration(long long n) override {
        if (n <= 3) return true;
        std::uniform_int_distribution<long long> dist(2, n - 2);
        long long a = dist(rng);
        return MathService::modPow(a, n - 1, n) == 1;
    }

public:
    FermatPrimalityTest() : rng(std::random_device{}()) {}
};
void runDemo() {
    FermatPrimalityTest fermatTest;
    long long numbers[] = { 2, 3, 5, 15, 97, 561 };
    double minProb = 0.99;

    std::cout << "=== Fermat Primality Test Demo ===" << std::endl;
    for (long long n : numbers) {
        bool isPrime = fermatTest.isPrime(n, minProb);
        std::cout << n << (isPrime ? " is probably prime" : " is composite") << std::endl;
    }
}
void runManualTests() {
    FermatPrimalityTest test;
    assert(test.isPrime(2, 0.9) == true);
    assert(test.isPrime(3, 0.9) == true);
    assert(test.isPrime(4, 0.9) == false);
    assert(test.isPrime(5, 0.9) == true);
    assert(test.isPrime(15, 0.9) == false);

    std::cout << "All manual Fermat tests passed!" << std::endl;
}

int main() {
    runDemo();
    runManualTests();
    return 0;
}

