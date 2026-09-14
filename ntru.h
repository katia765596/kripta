#ifndef NTRU_H
#define NTRU_H
#include "polynomial.h"
#include <string>
#include <utility>
struct ntru_parameters
{
    size_t N; big_int p; big_int q;};//степень полиномов, малый модуль и больщой (степень двойки или простое)
class ntru_encrypt
{
public:
    using polynomial_type = polynomial;
    struct public_key
    {
        ntru_parameters params;
        polynomial_type h;//h=p*f_q*g mod q
    };
    struct private_key
    {
        ntru_parameters params;
        polynomial_type f;//полином секретный
        polynomial_type f_p;//обратный к f по модулю p
    };
    struct ciphertext
    {
        polynomial_type e;
    };
    static std::pair<public_key, private_key> generate_keys(
        const ntru_parameters& params);
    static ciphertext encrypt(
        const public_key& pub,
        const polynomial_type& message);
    static polynomial_type decrypt(
        const private_key& priv,
        const ciphertext& cipher);
    static polynomial_type message_from_string(
        const std::string& text,
        const ntru_parameters& params);
    static std::string message_to_string(
        const polynomial_type& poly,
        const ntru_parameters& params);
};
#endif
