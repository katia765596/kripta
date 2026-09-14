#include "des.h"
#include <stdexcept>
namespace
{
    const int ip[64] = { 58,50,42,34,26,18,10,2,60,52,44,36,28,20,12,4,62,54,46,38,30,22,14,6,64,56,48,40,32,24,16,8,57,49,41,33,25,17,9,1,59,51,43,35,27,19,11,3,61,53,45,37,29,21,13,5,63,55,47,39,31,23,15,7 };
    const int fp[64] = { 40,8,48,16,56,24,64,32,39,7,47,15,55,23,63,31,38,6,46,14,54,22,62,30,37,5,45,13,53,21,61,29,36,4,44,12,52,20,60,28,35,3,43,11,51,19,59,27,34,2,42,10,50,18,58,26,33,1,41,9,49,17,57,25 };
    const int pc1[56] = { 57,49,41,33,25,17,9,1,58,50,42,34,26,18,10,2,59,51,43,35,27,19,11,3,60,52,44,36,63,55,47,39,31,23,15,7,62,54,46,38,30,22,14,6,61,53,45,37,29,21,13,5,28,20,12,4 };
    const int pc2[48] = { 14,17,11,24,1,5,3,28,15,6,21,10,23,19,12,4,26,8,16,7,27,20,13,2,41,52,31,37,47,55,30,40,51,45,33,48,44,49,39,56,34,53,46,42,50,36,29,32 };
    const int expansion[48] = { 32,1,2,3,4,5,4,5,6,7,8,9,8,9,10,11,12,13,12,13,14,15,16,17,16,17,18,19,20,21,20,21,22,23,24,25,24,25,26,27,28,29,28,29,30,31,32,1 };
    const int pbox[32] = { 16,7,20,21,29,12,28,17,1,15,23,26,5,18,31,10,2,8,24,14,32,27,3,9,19,13,30,6,22,11,4,25 };
    const int shifts[16] = { 1,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1 };
    const int sbox[8][64] = {
    {14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7,0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8,4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0,15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13},
    {15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10,3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5,0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15,13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9},
    {10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8,13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1,13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7,1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12},
    {7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15,13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9,10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4,3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14},
    {2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9,14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6,4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14,11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3},
    {12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11,10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8,9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6,4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13},
    {4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1,13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6,1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2,6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12},
    {13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7,1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2,7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8,2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11}
    };
}
uint64_t des::permute(uint64_t value, const int* table, int input_bits, int output_bits)
{
    uint64_t result = 0;
    for (int i = 0; i < output_bits; ++i)
        result = (result << 1) | ((value >> (input_bits - table[i])) & 1ULL);//сдвиг рез влево на 1 бит, освобождаем место, ищем смещение от старшего бита, извлекаем бит по указ смещению добавл в рез
    return result;
}
uint64_t des::load64(const byte_array& data)
{
    uint64_t value = 0;
    for (size_t i = 0; i < 8; ++i)
        value = (value << 8) | data[i];//первый байт старший
    return value;
}
void des::store64(byte_array& data, uint64_t value)
{
    for (int i = 7; i >= 0; --i)
    {
        data[i] = static_cast<uint8_t>(value & 0xff);
        value >>= 8;
    }
}
void des::generate_keys(const byte_array& key)//генерирует 16 раунд ключей на основе 8байтового
{
    uint64_t permuted = permute(load64(key), pc1, 64, 56);
    uint32_t left = static_cast<uint32_t>(permuted >> 28) & 0x0fffffff;//старшие 28 бит
    uint32_t right = static_cast<uint32_t>(permuted) & 0x0fffffff;
    for (int i = 0; i < 16; ++i)
    {
        left = ((left << shifts[i]) | (left >> (28 - shifts[i]))) & 0x0fffffff;//биты вышедшие за пределы переносятся вправо после сдвига влево
        right = ((right << shifts[i]) | (right >> (28 - shifts[i]))) & 0x0fffffff;
        round_keys[i] = permute((static_cast<uint64_t>(left) << 28) | right, pc2, 56, 48);//до 48 бит и сохр как раунд кей
    }
}
uint32_t des::feistel(uint32_t right, uint64_t key)
{
    uint64_t expanded = permute(right, expansion, 32, 48) ^ key;//расширяет правую половину до 48 бит и хорит с раунд ключом
    uint32_t result = 0;
    for (int i = 0; i < 8; ++i)//проход по 8 s-блокам
    {
        int six = static_cast<int>((expanded >> (42 - 6 * i)) & 0x3f);//излекает 6 бит
        int row = ((six & 0x20) >> 4) | (six & 1);//номер строки в s-блоке - старрший бит сдвигается на 4 (0 или 1 стал) и OR с младишм - дает номер строки (0..3)
        int column = (six >> 1) & 0x0f;//средние 4 бита (сдвиг вправо на 1)
        result = (result << 4) | sbox[i][row * 16 + column];//доб в рез 4-битное знач из s-блока (индекс:строка*16+столбец)
    }
    return static_cast<uint32_t>(permute(result, pbox, 32, 32));
}
byte_array des::crypt(const byte_array& block, bool decrypt) const
{
    if (block.size() != 8)
        throw std::invalid_argument("DES block size must be 8");
    uint64_t value = permute(load64(block), ip, 64, 64);//нач перестановка IP
    uint32_t left = static_cast<uint32_t>(value >> 32);
    uint32_t right = static_cast<uint32_t>(value);
    for (int i = 0; i < 16; ++i)
    {
        uint32_t next = left ^ feistel(right, round_keys[decrypt ? 15 - i : i]);//вычисляет новое знач для правой половины left XOR f(right, round_key), индекс раунд ключа зав от направления
        left = right;
        right = next;
    }
    uint64_t preoutput = (static_cast<uint64_t>(right) << 32) | left;
    byte_array result(8);
    store64(result, permute(preoutput, fp, 64, 64));//применяет конечную перестанову и сохр рез в байтовый массив
    return result;
}
des::des(const byte_array& key)
{
    if (key.size() != 8)
        throw std::invalid_argument("DES key size must be 8");
    generate_keys(key);
}
byte_array des::encrypt_block(const byte_array& block) const { return crypt(block, false); }
byte_array des::decrypt_block(const byte_array& block) const { return crypt(block, true); }
size_t des::get_block_size() const { return 8; }
size_t des::get_key_size() const { return 8; }
