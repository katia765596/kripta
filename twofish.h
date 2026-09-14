#ifndef TWOFISH_H
#define TWOFISH_H
#include "byte_array.h"
#include "block_cipher.h"
#include <cstddef>
#include <cstdint>
class twofish : public block_cipher
{
public:
    explicit twofish(const byte_array& key);
    byte_array encrypt_block(const byte_array& block) const override;
    byte_array decrypt_block(const byte_array& block) const override;
    size_t get_block_size() const override;
    size_t get_key_size() const;
    size_t get_round_count() const;
private:
    uint32_t round_keys[40];//40 раундовых ключей для 16 раундов
    uint32_t s_boxes[4][256];//4 s-блока по 256 элементов
    size_t key_size;
    static uint32_t rotl32(uint32_t value, unsigned int shift);
    static uint32_t rotr32(uint32_t value, unsigned int shift);
    static uint8_t gf_multiply(uint8_t a, uint8_t b, uint16_t polynomial);
    static uint8_t q_permutation(uint8_t value, bool q1);//q-перестановка (квадратичная подстановка) над байтом для s-блоков, q1 выбирает одну из двух перестановок q0 или q1
    static uint32_t load32(const byte_array& data, size_t offset);//загружает 32-битное из байтового массива начиная с указанноого смещения
    static void store32(byte_array& data, size_t offset, uint32_t value);
    static uint32_t h(uint32_t value, const uint32_t* key, size_t words);//реализует осн нелинейный преобраз
    static uint32_t h_full(uint32_t value, const uint32_t* key, size_t words);//полная версия h до шаг с испол маски и функции mds
    static uint32_t reed_solomon(uint32_t high, uint32_t low);//преобр Рида-Соломона при генерации s-блоков из ключа
    void generate_keys(const byte_array& key);
    uint32_t g(uint32_t value) const;
    uint32_t g0(uint32_t value) const;//испол в дешифровании для обратного преобраз
};
#endif
