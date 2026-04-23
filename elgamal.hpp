#ifndef ELGAMAL_HPP
#define ELGAMAL_HPP
#include "bigint.hpp"
#include <utility>
class ElGamal {
public:
    struct PublicKey { BigInt p, g, y; };
    struct PrivateKey { BigInt p, g, x; };
    struct Ciphertext { BigInt a, b; };
    static std::pair<PublicKey, PrivateKey> generate_keys(int bits);
    static Ciphertext encrypt(const PublicKey& pub, const BigInt& m);
    static BigInt decrypt(const PrivateKey& priv, const Ciphertext& cipher);
};
#endif 