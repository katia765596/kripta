#include "elgamal.hpp"
#include "ntru.hpp"
#include <iostream>
#include "elgamal.hpp"
#include "ntru.hpp"
#include <iostream>
#include <cassert>
void test_bigint() {
    using namespace std;
    BigInt a = 123456789;
    BigInt b("987654321");
    BigInt c = a + b;
    assert(c.to_string() == "1111111110");
    BigInt d = a * b;
    assert(d.to_string() == "121932631112635269");
    BigInt e = BigInt::mod_pow(BigInt(3), BigInt(10), BigInt(1000));
    assert(e == 49);
    BigInt p = BigInt::generate_prime(16);
    assert(BigInt::is_prime(p));
    cout << "[BigInt] tests passed\n";
}
void test_polynomial() {
    Polynomial a({ BigInt(1), BigInt(2), BigInt(3) });
    Polynomial b({ BigInt(0), BigInt(1), BigInt(1) });
    auto sum = a + b;
    assert(sum[0] == 1 && sum[1] == 3 && sum[2] == 4);
    auto prod = a * b;
    assert(prod.degree() == 4);
    auto mul_mod = a.mul_mod(b, BigInt(100), 3);
    assert(mul_mod[0] == 5 && mul_mod[1] == 4 && mul_mod[2] == 3);
    std::cout << "[Polynomial] tests passed\n";
}
void test_elgamal() {
    auto keys = ElGamal::generate_keys(32);
    BigInt msg = 12345;
    auto ct = ElGamal::encrypt(keys.first, msg);
    auto dec = ElGamal::decrypt(keys.second, ct);
    assert(msg == dec);
    std::cout << "[ElGamal] tests passed\n";
}
void test_ntru() {
    NTRUParameters params{ 11, BigInt(3), BigInt(61) };
    std::cout << "Generating NTRU keys..." << std::endl;
    auto keys = NTRUEncrypt::generate_keys(params);
    std::cout << "Keys generated." << std::endl;
    Polynomial msg({ BigInt(1), BigInt(2), BigInt(0), BigInt(1) });
    msg.resize(params.N);
    std::cout << "Original message polynomial: " << msg.to_string() << std::endl;
    auto ct = NTRUEncrypt::encrypt(keys.first, msg);
    std::cout << "Encrypted." << std::endl;
    auto dec = NTRUEncrypt::decrypt(keys.second, ct);
    std::cout << "Decrypted polynomial: " << dec.to_string() << std::endl;
    assert(msg == dec);
    std::cout << "Encryption/decryption roundtrip OK." << std::endl;
    std::string text = "Hello";
    auto mpoly = NTRUEncrypt::message_from_string(text, params);
    auto recovered = NTRUEncrypt::message_to_string(mpoly, params);
    assert(text.substr(0, params.N) == recovered);
    std::cout << "[NTRUEncrypt] tests passed\n";
}
int main() {
    test_bigint();
    test_polynomial();
    test_elgamal();
    test_ntru();
    std::cout << "\n ElGamal Demo \n";
    auto eg = ElGamal::generate_keys(64);
    BigInt msg = 0x123456789ABCDEF;
    auto ct = ElGamal::encrypt(eg.first, msg);
    auto dec = ElGamal::decrypt(eg.second, ct);
    std::cout << "Message: " << msg.to_string() << "\nDecrypted: " << dec.to_string() << "\n";
    std::cout << "\ NTRUEncrypt Demo \n";
    NTRUParameters ntru_params{ 11, BigInt(3), BigInt(61) };
    auto ntru = NTRUEncrypt::generate_keys(ntru_params);
    std::string text = "NTRU test";
    auto mpoly = NTRUEncrypt::message_from_string(text, ntru_params);
    auto ntru_ct = NTRUEncrypt::encrypt(ntru.first, mpoly);
    auto ntru_dec = NTRUEncrypt::decrypt(ntru.second, ntru_ct);
    std::string recovered = NTRUEncrypt::message_to_string(ntru_dec, ntru_params);
    std::cout << "Original: " << text << "\nRecovered: " << recovered << "\n";
    return 0;
}