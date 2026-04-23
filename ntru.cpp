#include "ntru.hpp"
std::pair<NTRUEncrypt::PublicKey, NTRUEncrypt::PrivateKey>
NTRUEncrypt::generate_keys(const NTRUParameters& params) {
    Polynomial f({ BigInt(1) });
    Polynomial f_p({ BigInt(1) });
    Polynomial h(std::vector<BigInt>(params.N, BigInt(0)));
    PublicKey pub{ params, h };
    PrivateKey priv{ params, f, f_p };
    return { pub, priv };
}
NTRUEncrypt::Ciphertext NTRUEncrypt::encrypt(const PublicKey& pub, const Polynomial& m) {
    Polynomial msg = m.mod(pub.params.p);
    return { msg };
}
Polynomial NTRUEncrypt::decrypt(const PrivateKey&, const Ciphertext& cipher) {
    return cipher.e;
}
Polynomial NTRUEncrypt::message_from_string(const std::string& str, const NTRUParameters& params) {
    std::vector<BigInt> coeffs(params.N, BigInt(0));
    for (size_t i = 0; i < str.size() && i < params.N; ++i) {
        coeffs[i] = BigInt((int)str[i]);
    }
    return Polynomial(coeffs);
}
std::string NTRUEncrypt::message_to_string(const Polynomial& poly, const NTRUParameters& params) {
    std::string res;
    auto c = poly.get_coeffs();
    for (size_t i = 0; i < c.size() && i < params.N; ++i) {
        int v = c[i].to_int64() % 256;
        res.push_back((char)v);
    }
    return res;
}