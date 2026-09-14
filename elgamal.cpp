#include "elgamal.h"
#include <stdexcept>
#include <vector>
namespace//анонимное пространство имен, все внутри имеют внутр связь (доступно в этом файле)
{
    std::vector<big_int> factorize(const big_int& n)//возв вектор простых множителей n для нахождения первообразного корня
    {
        std::vector<big_int> factors;
        big_int m = n;
        big_int d = 2;
        while (d * d <= m)//если делитель превысит корень из m, то оставшееся m будет простым
        {
            if (m % d == 0)
            {
                factors.push_back(d);
                while (m % d == 0)
                {
                    m /= d;
                }
            }
            if (d == 2)
            {
                d = 3;
            }
            else
            {
                d += 2;//только нечетные
            }
        }
        if (m > 1)
        {
            factors.push_back(m);//m посл простое 
        }
        return factors;
    }
}
std::pair<elgamal::public_key, elgamal::private_key>
elgamal::generate_keys(int bits)
{
    if (bits < 16)
    {
        throw std::invalid_argument("ElGamal requires at least 16 bits");
    }
    const big_int p = generate_prime(bits);
    const std::vector<big_int> factors = factorize(p - 1);//раскладываем на простые, это нужно для проверки, что кандидат в генераторы g явл первообразным корнем по модулю p
    big_int g = 2;
    while (true)
    {
        bool valid = true;
        for (const big_int& q : factors)
        {
            if (mod_pow(g, (p - 1) / q, p) == 1)//каждый дел удов g^((p-1)/q) mod p != 1
            {
                valid = false;
                break;
            }
        }
        if (valid)
        {
            break;
        }
        ++g;
    }
    const big_int x = random_bigint(2, p - 2);//1 и p-1 дают тривиал варианты при р-1 юзаем малую т ферма
    const big_int y = mod_pow(g, x, p);
    return {
        public_key{p, g, y},
        private_key{p, g, x}
    };
}
elgamal::ciphertext elgamal::encrypt(
    const public_key& pub,
    const big_int& message)
{
    if (message < 0 || message >= pub.p)
    {
        throw std::invalid_argument("message out of range");
    }
    const big_int k = random_bigint(2, pub.p - 2);//секретный ключ при шифровании каждого соо новый
    const big_int a = mod_pow(pub.g, k, pub.p);
    const big_int b = (message * mod_pow(pub.y, k, pub.p)) % pub.p;
    return { a, b };
}
big_int elgamal::decrypt(
    const private_key& priv,
    const ciphertext& cipher)
{
    if (cipher.a <= 0 || cipher.a >= priv.p ||
        cipher.b < 0 || cipher.b >= priv.p)
    {
        throw std::invalid_argument("invalid ciphertext");
    }
    const big_int shared_secret =
        mod_pow(cipher.a, priv.x, priv.p);
    return (cipher.b * mod_inverse(shared_secret, priv.p)) % priv.p;//назодим обратный для секрета по модулю p и умножаем на б и берем остаток
}
