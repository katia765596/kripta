#include "des.h"
#include "des_tables.h"
#include <cstring>
#include <stdexcept>
std::vector<std::vector<uint8_t>> des::des_key_schedule::expand_key(const std::vector<uint8_t>& key) {
    if (key.size() != 8) throw std::invalid_argument("des key must be 8 bytes");
    uint64_t k = bytes_to_uint64(key);
    uint64_t permuted = permute_64(k, des_pc1, 56); //перестановка PC=1 выбир 56из 64
    uint64_t c = (permuted >> 28) & 0x0FFFFFFF;//л половина (старшие биты 28-55) + маска
    uint64_t d = permuted & 0x0FFFFFFF;
    std::vector<std::vector<uint8_t>> round_keys;
    round_keys.reserve(16);
    for (int i = 0; i < 16; ++i) {
        int shift = des_shifts[i];//для каждого раунда свой сдвигз
        c = ((c << shift) | (c >> (28 - shift))) & 0x0FFFFFFF;//сдвиг влево и биты которые вышли за границу переносятся в младшие объединяем
        d = ((d << shift) | (d >> (28 - shift))) & 0x0FFFFFFF;
        uint64_t cd = (c << 28) | d;//объединяем ИЛИ побитовое в 56 бит (с в старшие 28)
        uint64_t k48 = permute_64(cd, des_pc2, 48);//PC=2 48 бит из 56 получаем раунд ключ
        round_keys.push_back(uint64_to_bytes(k48, 6));
    }
    return round_keys;
}
uint64_t des::des_key_schedule::bytes_to_uint64(const std::vector<uint8_t>& bytes) {
    uint64_t val = 0;
    for (size_t i = 0; i < 8; ++i) val = (val << 8) | bytes[i];//первый байт старший
    return val;
}
std::vector<uint8_t> des::des_key_schedule::uint64_to_bytes(uint64_t val, int len) {//64 битное знач и кол-во байт которые нужно извель
    std::vector<uint8_t> res(len);
    for (int i = len - 1; i >= 0; --i) {//берем байт через маску и записываем в рез (посл байт будет младший) и сдвиш вправо на 8 бит
        res[i] = val & 0xFF;
        val >>= 8;
    }
    return res;
}
uint64_t des::des_key_schedule::permute_64(uint64_t val, const int* table, int bits) {//из вал выбираем биты юзая массив где их порядок индексов в вых нач, битс- кол-во бит в вых знач
    uint64_t result = 0;
    for (int i = 0; i < bits; ++i) {
        int pos = table[i] - 1;//поз бита из табл в 0-индексации
        uint64_t bit = (val >> (63 - pos)) & 1;//бит в мл поз извл его (0 или 1)
        result = (result << 1) | bit;//сдвиг влеов (место для нового бита) и добавляем извлеченный с побитовым ИЛИ
    }
    return result;
}
std::vector<uint8_t> des::des_feistel_round::round_function(const std::vector<uint8_t>& block, const std::vector<uint8_t>& round_key) {
    if (block.size() != 4 || round_key.size() != 6) throw std::invalid_argument("des round input size error");
    uint32_t r = bytes_to_uint32(block);//п половина в 32 бита
    uint64_t r48 = expand_32_to_48(r);
    uint64_t k48 = bytes_to_uint48(round_key);//6 байт в 48 бит старшие нулевфе
    uint64_t x = r48 ^ k48;
    uint32_t output = 0;
    for (int i = 0; i < 8; ++i) {
        int six_bits = (x >> (42 - 6 * i)) & 0x3F;//извл 6 бит из х. сдвиг позволяет брать биты с поз 42 вниз для i=0 битв (47..42) и т д
        int row = ((six_bits & 0x20) >> 4) | (six_bits & 0x01);//номер строки в s-блоке 1 и последний бит 6 битного (сдвиг старшего на 4 поз,чтобы он стал вторым (0 или 2) затем или и получаем 0-3 номер строки
        int col = (six_bits >> 1) & 0x0F;//средние с 4го по 1-й - столбец, сдвиг на 1 убираем мл бит маска оставляет 4
        int sval = des_sbox[i][row * 16 + col];//вычисляем индекс (знач на пересечении) получаем 4-битное знач (0..15)
        output = (output << 4) | sval;
    }
    uint32_t permuted = permute_32(output, des_p_perm, 32);//P-перестановка 32 бита мешаем
    return uint32_to_bytes(permuted);//32 бита в 4 байта
}
uint32_t des::des_feistel_round::bytes_to_uint32(const std::vector<uint8_t>& bytes) {
    uint32_t val = 0;
    for (size_t i = 0; i < 4; ++i) val = (val << 8) | bytes[i];
    return val;
}
uint64_t des::des_feistel_round::bytes_to_uint48(const std::vector<uint8_t>& bytes) {
    uint64_t val = 0;
    for (size_t i = 0; i < 6; ++i) val = (val << 8) | bytes[i];
    return val;
}
std::vector<uint8_t> des::des_feistel_round::uint32_to_bytes(uint32_t val) {
    std::vector<uint8_t> res(4);
    for (int i = 3; i >= 0; --i) {
        res[i] = val & 0xFF;
        val >>= 8;
    }
    return res;
}
uint64_t des::des_feistel_round::expand_32_to_48(uint32_t val) {
    uint64_t result = 0;
    for (int i = 0; i < 48; ++i) {
        int pos = des_expansion[i] - 1;
        uint64_t bit = (val >> (31 - pos)) & 1;
        result = (result << 1) | bit;
    }
    return result;
}
uint32_t des::des_feistel_round::permute_32(uint32_t val, const int* table, int bits) {
    uint32_t result = 0;
    for (int i = 0; i < bits; ++i) {
        int pos = table[i] - 1;
        uint32_t bit = (val >> (31 - pos)) & 1;
        result = (result << 1) | bit;
    }
    return result;
}
des::des()
    : feistel_wrapper(new des_key_schedule(), new des_feistel_round(), 16),
    ks_holder(new des_key_schedule()),//констр иниц баз класс, констр хранит указ на динам объекты влож класс генер ключей и раунд функции
    fr_holder(new des_feistel_round()) {//баз класс отвечает за удал указ хранит как сырые указ
}
std::vector<uint8_t> des::pre_processing(const std::vector<uint8_t>& block) {//переопр вирт ф
    uint64_t b = bytes_to_uint64(block);
    uint64_t ip = permute_64(b, des_initial_perm, 64);//IP-перестановка перед 1 раундом
    return uint64_to_bytes(ip);
}
std::vector<uint8_t> des::post_processing(const std::vector<uint8_t>& block) {
    uint64_t b = bytes_to_uint64(block);
    uint64_t fp = permute_64(b, des_final_perm, 64);//FP-перестановка обр к ip после всех раундов
    return uint64_to_bytes(fp);
}
uint64_t des::bytes_to_uint64(const std::vector<uint8_t>& bytes) {
    uint64_t val = 0;
    for (size_t i = 0; i < 8; ++i) val = (val << 8) | bytes[i];
    return val;
}
std::vector<uint8_t> des::uint64_to_bytes(uint64_t val) {
    std::vector<uint8_t> res(8);
    for (int i = 7; i >= 0; --i) {
        res[i] = val & 0xFF;
        val >>= 8;
    }
    return res;
}
uint64_t des::permute_64(uint64_t val, const int* table, int bits) {
    uint64_t result = 0;
    for (int i = 0; i < bits; ++i) {
        int pos = table[i] - 1;
        uint64_t bit = (val >> (63 - pos)) & 1;
        result = (result << 1) | bit;
    }
    return result;
}