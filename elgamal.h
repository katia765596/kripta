#ifndef ELGAMAL_H
#define ELGAMAL_H
#include "big_int.h"
#include <utility>
class elgamal
{
public:
    struct public_key
    {
        big_int p;//простое чмсло
        big_int g;//первообразный корень по модулю p
        big_int y;//y=g^x mod p,x-секретный ключ
    };
    struct private_key
    {
        big_int p;
        big_int g;
        big_int x;
    };
    struct ciphertext
    {
        big_int a;//a=g^k mod p
        big_int b;//b=y^k*m mod p
    };
    static std::pair<public_key, private_key> generate_keys(int bits);
    static ciphertext encrypt(
        const public_key& pub,
        const big_int& message);
    static big_int decrypt(
        const private_key& priv,
        const ciphertext& cipher);
};
#endif
