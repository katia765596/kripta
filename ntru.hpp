#ifndef NTRU_HPP
#define NTRU_HPP
#include "polynomial.hpp"
#include <utility>
#include <string>
struct NTRUParameters {
    size_t N;
    BigInt p, q;
};
class NTRUEncrypt {
public:
    using Polynomial = ::Polynomial;
    struct PublicKey { NTRUParameters params; Polynomial h; };
    struct PrivateKey { NTRUParameters params; Polynomial f, f_p; };
    struct Ciphertext { Polynomial e; };
    static std::pair<PublicKey, PrivateKey> generate_keys(const NTRUParameters& params);
    static Ciphertext encrypt(const PublicKey& pub, const Polynomial& m);
    static Polynomial decrypt(const PrivateKey& priv, const Ciphertext& cipher);
    static Polynomial message_from_string(const std::string& str, const NTRUParameters& params);
    static std::string message_to_string(const Polynomial& poly, const NTRUParameters& params);
};
#endif 