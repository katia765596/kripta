#include "elgamal.h"
#include "ntru.h"
#include <cassert>
#include <iostream>
#include <string>
void test_bigint()
{
    const big_int a = 123456789;
    const big_int b("987654321");
    assert(a + b == big_int("1111111110"));
    assert(a * b == big_int("121932631112635269"));
    assert(mod_pow(3, 10, 1000) == 49);
    assert(is_prime(generate_prime(16)));
    std::cout << "big_int ok\n";
}
void test_polynomial()
{
    const polynomial a({ 1, 2, 3 });
    const polynomial b({ 0, 1, 1 });
    const polynomial sum = a + b;
    assert(sum[0] == 1);
    assert(sum[1] == 3);
    assert(sum[2] == 4);
    assert((a * b).degree() == 4);
    const polynomial product_mod =
        a.mul_mod(b, 100, 3);//3 - размерность N (степерь обрезания)
    assert(product_mod[0] == 5);
    assert(product_mod[1] == 4);
    assert(product_mod[2] == 3);
    std::cout << "polynomial ok\n";
}
void test_elgamal()
{
    const auto keys = elgamal::generate_keys(32);
    const big_int message = 12345;
    const auto cipher =
        elgamal::encrypt(keys.first, message);
    const big_int decrypted =
        elgamal::decrypt(keys.second, cipher);
    assert(message == decrypted);
    std::cout << "elgamal ok\n";
}
void test_ntru()
{
    const ntru_parameters params{ 11, 3, 12289 };
    std::cout << "generating NTRU keys...\n";
    const auto keys = ntru_encrypt::generate_keys(params);
    polynomial message({ 1, 2, 0, 1 });//полином сообщение 1+2x+0x^2+1x^3
    message.resize(params.N);//все полиномы должны иметь степень меньше N
    const auto cipher =
        ntru_encrypt::encrypt(keys.first, message);
    const auto decrypted =
        ntru_encrypt::decrypt(keys.second, cipher);
    assert(message == decrypted);
    std::cout << "ntru polynomial ok\n";
    const ntru_parameters text_params{ 11, 257, 65537 };
    const auto text_keys =
        ntru_encrypt::generate_keys(text_params);
    const std::string text = "Hello";
    const auto message_poly =
        ntru_encrypt::message_from_string(
            text,
            text_params);
    const auto text_cipher =
        ntru_encrypt::encrypt(
            text_keys.first,
            message_poly);
    const auto text_decrypted =
        ntru_encrypt::decrypt(
            text_keys.second,
            text_cipher);
    const std::string recovered =
        ntru_encrypt::message_to_string(
            text_decrypted,
            text_params);
    assert(recovered == text);
    std::cout << "ntru text ok: "
        << recovered << '\n';
}
int main()
{
    try
    {
        test_bigint();
        test_polynomial();
        test_elgamal();
        test_ntru();
        std::cout << "\ntests passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
