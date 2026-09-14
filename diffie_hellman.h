#ifndef DIFFIE_HELLMAN_H
#define DIFFIE_HELLMAN_H
#include <cstdint>
#include <cstddef>
#include "byte_array.h"
class diffie_hellman//протокол диффи хелмана для обмена ключами
{
public:
    diffie_hellman(uint64_t prime, uint64_t generator, uint64_t private_key);//констр принимает большое простое число(модуль),генератор (первообразный корень по модулю),секретный ключ участнника,вычисляет откр ключ ген^приват ключ мод прайм
    uint64_t get_prime() const;
    uint64_t get_generator() const;
    uint64_t get_private_key() const;
    uint64_t get_public_key() const;
    uint64_t make_shared_key(uint64_t other_public_key) const;//вычисляет общий секрет ключ на основе откр ключа другой стороны, откр^приват ключ мод прайм
    byte_array make_key(uint64_t other_public_key, size_t size) const;//вычисляет общий секрет и преобр в байтовый массив
private:
    uint64_t prime;
    uint64_t generator;
    uint64_t private_key;
    uint64_t public_key;
    static uint64_t multiply_mod(uint64_t a, uint64_t b, uint64_t mod);
    static uint64_t power_mod(uint64_t a, uint64_t e, uint64_t mod);
};
#endif
