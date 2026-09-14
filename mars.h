#ifndef MARS_H
#define MARS_H
#include "byte_array.h"
#include "block_cipher.h"//mars наследуется
#include <cstddef>
#include <cstdint>
class mars : public block_cipher
{
public:
    explicit mars(const byte_array& key);
    byte_array encrypt_block(const byte_array& block) const override;
    byte_array decrypt_block(const byte_array& block) const override;
    size_t get_block_size() const override;//размер блока фиксирован 16 байт
    size_t get_key_size() const;
    size_t get_round_count() const;//кол-во раундов зависит от размера ключа
private:
    uint32_t round_keys[40];//массив из 40 32-битных раунд ключей
    size_t key_size;
    static const uint32_t s_box[512];
    static uint32_t rotl32(uint32_t value, unsigned int shift);//цикл сдвиг 320битного значея влево
    static uint32_t rotr32(uint32_t value, unsigned int shift);
    static uint32_t load32(const byte_array& data, size_t offset);//загружает 32-битное знач из массива байтов по смещению
    static void store32(byte_array& data, size_t offset, uint32_t value);//сохрает 32 битное знач в массив байтов по смещению
    void generate_keys(const byte_array& key);
};
#endif
