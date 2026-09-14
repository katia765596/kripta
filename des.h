#ifndef DES_H
#define DES_H
#include "byte_array.h"
#include "block_cipher.h"
#include <cstddef>
#include <cstdint>
class des : public block_cipher//размер блока 64 бита и ключ 56 бит
{
public:
    explicit des(const byte_array& key);
    byte_array encrypt_block(const byte_array& block) const override;
    byte_array decrypt_block(const byte_array& block) const override;
    size_t get_block_size() const override;
    size_t get_key_size() const;
private:
    uint64_t round_keys[16];//фактически 48 значащих бит
    static uint64_t permute(uint64_t value, const int* table, int input_bits, int output_bits);
    static uint32_t feistel(uint32_t right, uint64_t key);//раун функция фейстеля принимает правую половину 32 бита и раун ключ
    static uint64_t load64(const byte_array& data);//загружает 64битное из байтового массива в порядке бигэнджен
    static void store64(byte_array& data, uint64_t value);//сохр 64битное знач в байтовый массив
    void generate_keys(const byte_array& key);
    byte_array crypt(const byte_array& block, bool decrypt) const;//испол раунд ключи в прямом или обратное порядке
};
#endif
