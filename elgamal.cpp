#include "elgamal.hpp"
#include <vector>
static std::vector<BigInt> factorize(const BigInt& n) {
    std::vector<BigInt> factors;
    BigInt m = n;
    BigInt d = 2;
    while (d * d <= m) {
        if (m % d == 0) {
            factors.push_back(d);
            while (m % d == 0) m = m / d;
        }
        if (d == 2) d = 3;
        else d += 2;
    }
    if (m > 1) factors.push_back(m);
    return factors;
}
std::pair<ElGamal::PublicKey, ElGamal::PrivateKey> ElGamal::generate_keys(int bits) {
    BigInt p = BigInt::generate_prime(bits);
    BigInt phi = p - 1;
    auto prime_factors = factorize(phi);
    BigInt g = 2;
    while (true) {
        bool ok = true;
        for (const auto& q : prime_factors) {
            if (BigInt::mod_pow(g, phi / q, p) == 1) { ok = false; break; }
        }
        if (ok) break;
        ++g;
    }
    BigInt x = BigInt::random(2, p - 2);
    BigInt y = BigInt::mod_pow(g, x, p);
    return { {p, g, y}, {p, g, x} };
}
ElGamal::Ciphertext ElGamal::encrypt(const PublicKey& pub, const BigInt& m) {
    BigInt k = BigInt::random(2, pub.p - 2);
    BigInt a = BigInt::mod_pow(pub.g, k, pub.p);
    BigInt b = (m * BigInt::mod_pow(pub.y, k, pub.p)) % pub.p;
    return { a, b };
}
BigInt ElGamal::decrypt(const PrivateKey& priv, const Ciphertext& cipher) {
    BigInt s = BigInt::mod_pow(cipher.a, priv.x, priv.p);
    BigInt inv_s = BigInt::mod_inverse(s, priv.p);
    return (cipher.b * inv_s) % priv.p;
}