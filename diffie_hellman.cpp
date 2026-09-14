#include "diffie_hellman.h"
#include <stdexcept>
diffie_hellman::diffie_hellman(uint64_t prime_value, uint64_t generator_value, uint64_t private_key_value)
    : prime(prime_value), generator(generator_value), private_key(private_key_value), public_key(0)
{
    if (prime < 3 || generator <= 1 || generator >= prime || private_key == 0)//если приват=0 то открытый 1, тривиально
        throw std::invalid_argument("invalid Diffie-Hellman parameters");
    public_key = power_mod(generator, private_key, prime);
}
uint64_t diffie_hellman::multiply_mod(uint64_t a, uint64_t b, uint64_t mod)
{
    uint64_t result = 0;
    a %= mod;
    while (b != 0)
    {
        if (b & 1)
            result = result >= mod - a ? result - (mod - a) : result + a;//если младший бит установлен 1 то прибавляем а,если рез больше мод-а, то рез вычитает mod-a
        a = a >= mod - a ? a - (mod - a) : a + a;//удваивает а по мод
        b >>= 1;//переход к след биту
    }
    return result;
}
uint64_t diffie_hellman::power_mod(uint64_t a, uint64_t e, uint64_t mod)
{
    uint64_t result = 1 % mod;
    a %= mod;
    while (e != 0)
    {
        if (e & 1)
            result = multiply_mod(result, a, mod);
        a = multiply_mod(a, a, mod);
        e >>= 1;
    }
    return result;
}//геттеры методы доступа к приватным полям
uint64_t diffie_hellman::get_prime() const { return prime; }
uint64_t diffie_hellman::get_generator() const { return generator; }
uint64_t diffie_hellman::get_private_key() const { return private_key; }
uint64_t diffie_hellman::get_public_key() const { return public_key; }
uint64_t diffie_hellman::make_shared_key(uint64_t other_public_key) const
{
    if (other_public_key == 0 || other_public_key >= prime)
        throw std::invalid_argument("invalid public key");
    return power_mod(other_public_key, private_key, prime);
}//вычисляет общий секретный ключ по открытому ключу другой стороны
byte_array diffie_hellman::make_key(uint64_t other_public_key, size_t size) const
{//общий секрет в байт массив
    if (size == 0)
        return {};
    uint64_t shared = make_shared_key(other_public_key);
    byte_array result(size);
    for (size_t i = 0; i < size; ++i)
    {
        uint64_t value = power_mod(shared + static_cast<uint64_t>(i + 1), static_cast<uint64_t>(i + 3), prime);//если запросить сайз больше 8, то байтов не хватит, нужно генер больше псевдослуч,для каждого i дает новое псевдослуч 64битное
        result[i] = static_cast<uint8_t>((value ^ (shared >> ((i % 8) * 8))) & 0xff);//извлекает из shared байт с номером i%8, хор с value добавляет вклад исходного секрета
    }
    return result;
}
