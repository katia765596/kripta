#include "twofish.h"//симм блочный шифр с размером блока 128 бит и длиной ключа до 256 бит
#include <stdexcept>
namespace
{
    const uint8_t q_t[2][4][16] = {//два Q блока по 4 таблицы по 16 эл, с помощью табл преобр 8битные знач, основа для s-боксов
        {
            {8,1,7,13,6,15,3,2,0,11,5,9,14,12,10,4},
            {14,12,11,8,1,2,3,5,15,4,10,6,7,0,9,13},
            {11,10,5,14,6,13,9,0,12,8,15,3,2,4,7,1},
            {13,7,15,4,1,2,6,14,9,11,3,0,8,5,12,10}
        },
        {
            {2,8,11,13,15,7,6,14,3,1,9,4,0,10,12,5},
            {1,14,2,11,4,12,3,7,6,13,10,5,15,9,0,8},
            {4,12,7,5,1,6,9,10,0,14,13,8,2,11,3,15},
            {11,9,5,1,12,3,13,14,6,4,7,15,2,0,8,10}
        }
    };
    const uint8_t rs_matrix[4][8] = {//испол при генерации s-боксов из ключа
        {0x01,0xa4,0x55,0x87,0x5a,0x58,0xdb,0x9e},
        {0xa4,0x56,0x82,0xf3,0x1e,0xc6,0x68,0xe5},
        {0x02,0xa1,0xfc,0xc1,0x47,0xae,0x3d,0x19},
        {0xa4,0x55,0x87,0x5a,0x58,0xdb,0x9e,0x03}
    };
    const uint8_t mds_matrix[4][4] = {//испол в преобр h (часть шифр преобр) для обеспечения диффузии, коэф поле галуа
        {0x01,0xef,0x5b,0x5b},
        {0x5b,0xef,0xef,0x01},
        {0xef,0x5b,0x01,0xef},
        {0xef,0x01,0xef,0x5b}
    };
    uint8_t ror4(uint8_t value)
    {
        return static_cast<uint8_t>((value >> 1) | ((value & 1) << 3));//цикл сдвиг 8-битного знач (младшие 4 бита) вправо на 1 бит,извекаем младший бит исх знач, сдвигает мл бит на 3 поз влево, помещая его в бит 3 четвертый бит младших битов,объед
    }
    uint8_t gf_mul(uint8_t a, uint8_t b, uint16_t polynomial)//умн в поле
    {
        uint16_t result = 0;
        uint16_t value = b;
        while (a != 0)
        {
            if (a & 1)//если младший бит а равен 1,хорит рез с тек вал
                result ^= value;
            a >>= 1;//для обработки след бита
            value <<= 1;//влево на 1 бит умн на х
            if (value & 0x100)
                value ^= polynomial;//по модулю
        }
        return static_cast<uint8_t>(result);
    }
    uint8_t q_permutation_value(uint8_t value, bool q1)
    {
        uint8_t a0 = value >> 4;//4 бита входного байта старшие
        uint8_t b0 = value & 15;//младшие 4 бита
        uint8_t a1 = a0 ^ b0;
        uint8_t b1 = a0 ^ ror4(b0) ^ static_cast<uint8_t>((8 * a0) & 15);//хор а0, цикл сдвинутой вправо b0 и 
        uint8_t a2 = q_t[q1 ? 1 : 0][0][a1];//применяет первую табл индекс 0 к a1
        uint8_t b2 = q_t[q1 ? 1 : 0][1][b1];//применяет вторую табл к b1
        uint8_t a3 = a2 ^ b2;
        uint8_t b3 = a2 ^ ror4(b2) ^ static_cast<uint8_t>((8 * a2) & 15);//аналогично b1
        uint8_t a4 = q_t[q1 ? 1 : 0][2][a3];//третья
        uint8_t b4 = q_t[q1 ? 1 : 0][3][b3];//четвертая
        return static_cast<uint8_t>((b4 << 4) | a4);//возвр байт старшие и младшие
    }
    uint8_t h_byte_value(int k, int i, uint8_t x, uint8_t l0, uint8_t l1, uint8_t l2, uint8_t l3)
    {//h-функция для генер ключевых зав, выпол неск послед Q-перестановок с XOR с байтами ключ слов
        const bool q_stage[4][5] = {//матрица определ какую q-табл испол на каждом этапе для разных индексов и размера ключа
            {true, true, false, false, true},
            {false, true, true, false, false},
            {false, false, false, true, true},
            {true, false, true, true, false}
        };
        uint8_t l[4] = { l0, l1, l2, l3 };//4 байта ключ слов
        uint8_t v;
        if (k == 2)//ключ из двух слов
        {
            v = q_permutation_value(x, q_stage[i][2]);//Q-перестановка к х с выбором табл
        }
        else if (k == 3)
        {
            v = q_permutation_value(x, q_stage[i][1]);
            v = q_permutation_value(static_cast<uint8_t>(l[2] ^ v), q_stage[i][2]);//там еще хор с v
        }
        else
        {
            v = q_permutation_value(x, q_stage[i][0]);
            v = q_permutation_value(static_cast<uint8_t>(l[3] ^ v), q_stage[i][1]);
            v = q_permutation_value(static_cast<uint8_t>(l[2] ^ v), q_stage[i][2]);
        }
        v = q_permutation_value(static_cast<uint8_t>(l[1] ^ v), q_stage[i][3]);
        v = q_permutation_value(static_cast<uint8_t>(l[0] ^ v), q_stage[i][4]);
        return v;
    }
    uint32_t h_byte_mds(int k, int i, uint8_t x, uint8_t l0, uint8_t l1, uint8_t l2, uint8_t l3)
    {//применяет h-байтовую ф, а затем умножает рез на столбец MDS-матрицы (матр макс расстояния разделения) для получения 32 битного
        uint8_t v = h_byte_value(k, i, x, l0, l1, l2, l3);
        uint32_t result = 0;
        for (int row = 0; row < 4; ++row)
            result |= static_cast<uint32_t>(gf_mul(mds_matrix[row][i], v, 0x69)) << (8 * row);
        return result;//умн эл mds-матр на v в поле галуа с полиномм 0x69, рез каждого умн сдвигается влево на 8*row и накапл в рез
    }
    uint32_t h_word(int k, uint8_t x, uint32_t l0, uint32_t l1, uint32_t l2, uint32_t l3)
    {//Н-функция для 32бит слва x 
        uint8_t l[4][4] = {};//4 ключевые слова
        uint32_t words[4] = { l0,l1,l2,l3 };//массив из четырех 32битных слов ключа
        for (int w = 0; w < 4; ++w)
            for (int i = 0; i < 4; ++i)
                l[w][i] = static_cast<uint8_t>(words[w] >> (8 * i));//преоб каждое слово в 4 байта
        uint32_t result = 0;
        for (int i = 0; i < 4; ++i)
        {
            uint8_t v = h_byte_value(k, i, x, l[0][i], l[1][i], l[2][i], l[3][i]);
            for (int row = 0; row < 4; ++row)
                result ^= static_cast<uint32_t>(gf_mul(mds_matrix[row][i], v, 0x69)) << (8 * row);
        }
        return result;
    }
    uint32_t compute_s(uint32_t m1, uint32_t m2)//принимает два 32битных возвр 32битное, реал умн на матрицу рида для генерации s-блоков
    {
        uint32_t result = 0;
        for (int i = 0; i < 4; ++i)
        {
            uint8_t value =
                gf_mul(static_cast<uint8_t>(m1), rs_matrix[i][0], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m1 >> 8), rs_matrix[i][1], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m1 >> 16), rs_matrix[i][2], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m1 >> 24), rs_matrix[i][3], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m2), rs_matrix[i][4], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m2 >> 8), rs_matrix[i][5], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m2 >> 16), rs_matrix[i][6], 0x4d) ^
                gf_mul(static_cast<uint8_t>(m2 >> 24), rs_matrix[i][7], 0x4d);
            result |= static_cast<uint32_t>(value) << (8 * i);
        }
        return result;
    }//для каэжой строки i и каждого из 8 байтов выполн умн в поле галуа, на соотв эл rs-матр с полин 0x4D , рез хорится
}//возвр 32-бинтео составленное из 4 байт
uint32_t twofish::rotl32(uint32_t value, unsigned int shift)
{
    shift &= 31;
    return shift == 0 ? value : (value << shift) | (value >> (32 - shift));
}//сдвиг 32битного на шифт поз
uint32_t twofish::rotr32(uint32_t value, unsigned int shift)
{
    shift &= 31;
    return shift == 0 ? value : (value >> shift) | (value << (32 - shift));
}//сдвиг вправо
uint8_t twofish::gf_multiply(uint8_t a, uint8_t b, uint16_t polynomial)
{
    return gf_mul(a, b, polynomial);
}
uint8_t twofish::q_permutation(uint8_t value, bool q1)
{
    return q_permutation_value(value, q1);
}
uint32_t twofish::load32(const byte_array& data, size_t offset)
{//загружает 32 битное из массива по смещению
    return static_cast<uint32_t>(data[offset]) |
        (static_cast<uint32_t>(data[offset + 1]) << 8) |
        (static_cast<uint32_t>(data[offset + 2]) << 16) |
        (static_cast<uint32_t>(data[offset + 3]) << 24);
}
void twofish::store32(byte_array& data, size_t offset, uint32_t value)
{//сохраняет 32 битное в массив по смещению
    data[offset] = static_cast<uint8_t>(value);
    data[offset + 1] = static_cast<uint8_t>(value >> 8);
    data[offset + 2] = static_cast<uint8_t>(value >> 16);
    data[offset + 3] = static_cast<uint8_t>(value >> 24);
}
uint32_t twofish::reed_solomon(uint32_t high, uint32_t low)
{
    return compute_s(high, low);//генер s-боксов
}
uint32_t twofish::h_full(uint32_t value, const uint32_t* key, size_t words)
{
    uint32_t l0 = key[0];
    uint32_t l1 = key[1];
    uint32_t l2 = words > 2 ? key[2] : 0;
    uint32_t l3 = words > 3 ? key[3] : 0;
    return h_word(static_cast<int>(words), static_cast<uint8_t>(value), l0, l1, l2, l3);
}//вызывает h_wprd, передавая слва как инт млад байт, возвр 32битны 
uint32_t twofish::h(uint32_t value, const uint32_t* key, size_t words)
{
    return h_full(value, key, words);
}//обертка
void twofish::generate_keys(const byte_array& key)
{
    uint32_t m[8] = {};//массив из 8 32битных слво
    for (size_t i = 0; i < key.size() / 4; ++i)
        m[i] = load32(key, i * 4);//загружает ключ в массив m
    int k = static_cast<int>(key.size() / 8);//кол-во 64битных половин ключа
    for (int i = 0; i < 20; ++i)
    {//ген 40 раунд ключей по 2 на каждый из 20 раундов
        uint32_t a = h_word(k, static_cast<uint8_t>(2 * i), m[0], m[2], m[4], m[6]);
        uint32_t b = h_word(k, static_cast<uint8_t>(2 * i + 1), m[1], m[3], m[5], m[7]);
        b = rotl32(b, 8);//сдвиг б влево на 8 бит
        round_keys[2 * i] = a + b;//первый раунд ключ  сложение по модулю 2^32)
        round_keys[2 * i + 1] = rotl32(a + b + b, 9);//второй раунд ключ (a + 2*b) с циклическим сдвигом влево на 9 бит
    }
    uint32_t s[4] = {};//массив из 4 32битных для s-боксов
    for (int i = 0; i < k; ++i)
        s[k - 1 - i] = compute_s(m[2 * i], m[2 * i + 1]);//вычисляет s-слова в обр порядке для каждой пары ключ слов
    for (int i = 0; i < 4; ++i)//по sбоксам
        for (int j = 0; j < 256; ++j)
            s_boxes[i][j] = h_byte_mds(k, i, static_cast<uint8_t>(j), static_cast<uint8_t>(s[0] >> (8 * i)), static_cast<uint8_t>(s[1] >> (8 * i)), static_cast<uint8_t>(s[2] >> (8 * i)), static_cast<uint8_t>(s[3] >> (8 * i)));
}//вычисляет знач s-бокса, из каждого s-слова извлекается i бит и передается в h_byte_mds 
twofish::twofish(const byte_array& key)
    : key_size(key.size())
{
    if (key_size != 16 && key_size != 24 && key_size != 32)
        throw std::invalid_argument("invalid Twofish key size");
    generate_keys(key);//ген раунд ключей и s-боксов на основе переданного ключа
}
uint32_t twofish::g(uint32_t value) const
{
    return s_boxes[1][static_cast<uint8_t>(value)] ^//извлекает младший байт из value
        s_boxes[2][static_cast<uint8_t>(value >> 8)] ^
        s_boxes[3][static_cast<uint8_t>(value >> 16)] ^
        s_boxes[0][static_cast<uint8_t>(value >> 24)];
}//замена каждого байта через s-боксы с последующим хоr
uint32_t twofish::g0(uint32_t value) const
{
    return s_boxes[0][static_cast<uint8_t>(value)] ^
        s_boxes[1][static_cast<uint8_t>(value >> 8)] ^
        s_boxes[2][static_cast<uint8_t>(value >> 16)] ^
        s_boxes[3][static_cast<uint8_t>(value >> 24)];
}
byte_array twofish::encrypt_block(const byte_array& block) const
{
    if (block.size() != 16)
        throw std::invalid_argument("invalid block size");
    uint32_t r0 = load32(block, 0) ^ round_keys[0];//загрузка 4 32 битных слов из блока
    uint32_t r1 = load32(block, 4) ^ round_keys[1];
    uint32_t r2 = load32(block, 8) ^ round_keys[2];
    uint32_t r3 = load32(block, 12) ^ round_keys[3];
    for (int i = 0; i < 8; ++i)//осн цикл из 8 итер (16 раундов, но в коде сгруппированы по 2 подраунда на каждую итер)
    {
        uint32_t t1 = g(r1);//сдвиг влево на 8 бит
        uint32_t t0 = g0(r0) + t1;
        r3 = (t1 + t0 + round_keys[4 * i + 9]) ^ rotl32(r3, 1);
        r2 = rotr32((t0 + round_keys[4 * i + 8]) ^ r2, 1);
        t1 = g(r3);
        t0 = g0(r2) + t1;
        r1 = (t1 + t0 + round_keys[4 * i + 11]) ^ rotl32(r1, 1);
        r0 = rotr32((t0 + round_keys[4 * i + 10]) ^ r0, 1);
    }
    byte_array result(16);
    store32(result, 0, r2 ^ round_keys[4]);
    store32(result, 4, r3 ^ round_keys[5]);
    store32(result, 8, r0 ^ round_keys[6]);
    store32(result, 12, r1 ^ round_keys[7]);
    return result;
}
byte_array twofish::decrypt_block(const byte_array& block) const
{
    if (block.size() != 16)
        throw std::invalid_argument("invalid block size");
    uint32_t r0 = load32(block, 8) ^ round_keys[6];
    uint32_t r1 = load32(block, 12) ^ round_keys[7];
    uint32_t r2 = load32(block, 0) ^ round_keys[4];
    uint32_t r3 = load32(block, 4) ^ round_keys[5];
    for (int i = 0; i < 8; ++i)
    {
        uint32_t t1 = g(r3);
        uint32_t t0 = g0(r2) + t1;
        r1 = rotr32((t1 + t0 + round_keys[39 - 4 * i]) ^ r1, 1);
        r0 = (t0 + round_keys[38 - 4 * i]) ^ rotl32(r0, 1);
        t1 = g(r1);
        t0 = g0(r0) + t1;
        r3 = rotr32((t1 + t0 + round_keys[37 - 4 * i]) ^ r3, 1);
        r2 = (t0 + round_keys[36 - 4 * i]) ^ rotl32(r2, 1);
    }
    byte_array result(16);
    store32(result, 0, r0 ^ round_keys[0]);
    store32(result, 4, r1 ^ round_keys[1]);
    store32(result, 8, r2 ^ round_keys[2]);
    store32(result, 12, r3 ^ round_keys[3]);
    return result;
}
size_t twofish::get_block_size() const
{
    return 16;
}
size_t twofish::get_key_size() const
{
    return key_size;
}
size_t twofish::get_round_count() const
{
    return 16;
}
