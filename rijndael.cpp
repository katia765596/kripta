#include "rijndael.h"
#include <stdexcept>
uint8_t rijndael::multiply(uint8_t first, uint8_t second)//неприводимый x^8 + x^4 + x^3 + x + 1
{
    uint8_t result = 0;
    for (int i = 0; i < 8; ++i)//выполняем 8 раз по кол-ву бит в байте, на каждой итер обрабатывается один бит секонд
    {
        if (second & 1)//если младший бит 1 у секрнд
        {
            result ^= first;
        }
        bool high = (first & 0x80) != 0;//знач старшего бита 7-го переменной ферст перед сдвигом
        first <<= 1;//сдвиг ферст на 1бит влево (умн на х)
        if (high)//если бит до сдвига был равен 1, то после член степени 8, редуцируем по неприв полиному без старшего бита
        {
            first ^= 0x1b;
        }
        second >>= 1;//сдвиг вправо на 1 бит, чтобы на след итер проверить след бит
    }
    return result;
}
uint8_t rijndael::gf_inverse(uint8_t value)
{
    if (value == 0)
    {
        return 0;//в поле галуа обр эл для 0 нет, но в aes ноль
    }
    uint8_t result = 1;
    uint8_t base = value;
    unsigned int exponent = 254;//value^255 = 1
    while (exponent != 0)
    {
        if (exponent & 1)//мл бит 1
        {
            result = multiply(result, base);
        }
        base = multiply(base, base);
        exponent >>= 1;//вправо на 1 бит
    }
    return result;
}
uint8_t rijndael::rotl8(uint8_t value, unsigned int shift)
{
    shift %= 8;//огр шифт знач от 0 до 7 тк сдвиг на 8 равносилен сдвигу на остаток от деления на 8
    if (shift == 0)
    {
        return value;
    }
    return static_cast<uint8_t>(//часть выходящая за левую гран переносится в правую
        (value << shift) |//младшие заполняются 0
        (value >> (8 - shift))//старшие в конец
        );
}
rijndael::rijndael(
    size_t block_size_value,
    const byte_array& key_value
)
    : block_size(block_size_value),
    key_size(key_value.size()),
    nb(block_size_value / 4),
    nk(key_value.size() / 4),
    nr(0),
    key(key_value)
{
    if (block_size != 16 && block_size != 24 && block_size != 32)
    {
        throw std::invalid_argument("invalid Rijndael block size");
    }
    if (key_size != 16 && key_size != 24 && key_size != 32)
    {
        throw std::invalid_argument("invalid Rijndael key size");
    }
    nr = (nb > nk ? nb : nk) + 6;//кол-во рацндов макс(nb,nk)+6
    generate_s_boxes();
    generate_round_keys();
}
void rijndael::generate_s_boxes()
{
    s_box.resize(256);
    inv_s_box.resize(256);
    for (int i = 0; i < 256; ++i)
    {
        uint8_t value = static_cast<uint8_t>(i);//тек знач байта
        uint8_t inverse = gf_inverse(value);
        s_box[i] = static_cast<uint8_t>(
            inverse ^
            rotl8(inverse, 1) ^
            rotl8(inverse, 2) ^
            rotl8(inverse, 3) ^
            rotl8(inverse, 4) ^
            0x63
            );//афин преоб над инверс построение s-box
    }
    for (int i = 0; i < 256; ++i)//обр s-box 
    {
        inv_s_box[s_box[i]] = static_cast<uint8_t>(i);//для каждого i берем знач s_box[i] и по этому знач индексируем обр s-box,записывая туда i
    }
}
uint8_t rijndael::get_round_constant(size_t index) const//раунд константа для заданного индекса
{
    uint8_t result = 1;
    if (index == 0)
    {
        return 0;
    }
    for (size_t i = 1; i < index; ++i)
    {
        result = multiply(result, 2);
    }
    return result;//2^(index-1)
}
void rijndael::generate_round_keys()
{
    size_t words = nb * (nr + 1);//общ кол-во 4байтовых слов, необх для хранения всех раунд ключей еще 1 для начального
    round_keys.resize(words * 4);
    for (size_t i = 0; i < key_size; ++i)
    {
        round_keys[i] = key[i];
    }//первый раунд ключ раунд 0
    size_t generated = nk;//кол-во уже сгенерированных 4-байтовых слов изначально nk слов - исх ключ
    while (generated < words)
    {
        uint8_t temp[4];//4 байта, для хранения промежуточного слова
        size_t previous = (generated - 1) * 4;//смещение в байтах для посл сген слова
        for (size_t i = 0; i < 4; ++i)
        {
            temp[i] = round_keys[previous + i];
        }//копия 4 байт предыд слова в массив temp
        if (generated % nk == 0)
        {
            uint8_t value = temp[0];//сохр первый байт
            temp[0] = s_box[temp[1]];//преоб RotWord свдиг влево на 1 байт и замена байта через s-box
            temp[1] = s_box[temp[2]];
            temp[2] = s_box[temp[3]];
            temp[3] = s_box[value];
            temp[0] ^= get_round_constant(generated / nk);//хорим первый байт с раунд конст, соответ номеру тек блока
        }
        else if (nk > 6 && generated % nk == 4)//для 256битного
        {
            temp[0] = s_box[temp[0]];
            temp[1] = s_box[temp[1]];
            temp[2] = s_box[temp[2]];
            temp[3] = s_box[temp[3]];
        }
        size_t base = (generated - nk) * 4;
        size_t position = generated * 4;//куда записано новое слово
        for (size_t i = 0; i < 4; ++i)
        {
            round_keys[position + i] =
                round_keys[base + i] ^ temp[i];
        }//вычисляет новое слово как хор слова, находящегося на nk шагов назад и преобр тем
        ++generated;
    }
}
void rijndael::sub_bytes(byte_array& state) const
{
    for (size_t i = 0; i < state.size(); ++i)
    {
        state[i] = s_box[state[i]];//для каждого байта состояния заменяет его на соотв знач из S-бокс
    }
}
void rijndael::inv_sub_bytes(byte_array& state) const
{
    for (size_t i = 0; i < state.size(); ++i)
    {
        state[i] = inv_s_box[state[i]];
    }
}
void rijndael::shift_rows(byte_array& state) const
{
    byte_array result = state;//чтобы не перезаписывать исх состояние до завершения перестановки
    for (size_t row = 0; row < 4; ++row)
    {
        for (size_t column = 0; column < nb; ++column)
        {
            size_t source = (column + row) % nb;//индекс исх столбца для тек столбца column с учетом сдвига на row позиций влеов
            result[4 * column + row] = state[4 * source + row];//записть в рез значения из исх состояния, бай в столбце с и строке r находится по индексу 4*с+r
        }
    }
    state = result;
}
void rijndael::inv_shift_rows(byte_array& state) const
{
    byte_array result = state;
    for (size_t row = 0; row < 4; ++row)
    {
        for (size_t column = 0; column < nb; ++column)
        {
            size_t source = (column + nb - row) % nb;//для обр сдвига столбец вычисл -row эквиваленто сдвигу вправо на row поз
            result[4 * column + row] = state[4 * source + row];
        }
    }
    state = result;
}
void rijndael::mix_columns(byte_array& state) const
{//умн каждого столбца сост на фиксир матрицу
    for (size_t column = 0; column < nb; ++column)
    {
        size_t index = column * 4;//индекс в байтовом массиве, соответ началу тек столбца
        uint8_t a0 = state[index];
        uint8_t a1 = state[index + 1];
        uint8_t a2 = state[index + 2];
        uint8_t a3 = state[index + 3];
        state[index] =
            multiply(a0, 2) ^ multiply(a1, 3) ^ a2 ^ a3;
        state[index + 1] =
            a0 ^ multiply(a1, 2) ^ multiply(a2, 3) ^ a3;
        state[index + 2] =
            a0 ^ a1 ^ multiply(a2, 2) ^ multiply(a3, 3);
        state[index + 3] =
            multiply(a0, 3) ^ a1 ^ a2 ^ multiply(a3, 2);
    }
}
void rijndael::inv_mix_columns(byte_array& state) const
{
    for (size_t column = 0; column < nb; ++column)
    {
        size_t index = column * 4;
        uint8_t a0 = state[index];
        uint8_t a1 = state[index + 1];
        uint8_t a2 = state[index + 2];
        uint8_t a3 = state[index + 3];
        state[index] =
            multiply(a0, 14) ^ multiply(a1, 11) ^
            multiply(a2, 13) ^ multiply(a3, 9);
        state[index + 1] =
            multiply(a0, 9) ^ multiply(a1, 14) ^
            multiply(a2, 11) ^ multiply(a3, 13);
        state[index + 2] =
            multiply(a0, 13) ^ multiply(a1, 9) ^
            multiply(a2, 14) ^ multiply(a3, 11);
        state[index + 3] =
            multiply(a0, 11) ^ multiply(a1, 13) ^
            multiply(a2, 9) ^ multiply(a3, 14);
    }
}//хор-сложение в поле гф(2)
void rijndael::add_round_key(byte_array& state, size_t round) const
{//хор состояния с раунд ключом для заданоного раудна
    size_t offset = round * block_size;//вычисляет смещение в массиве раунд ключей, где начинается ключ для данного раунда
    for (size_t i = 0; i < block_size; ++i)
    {
        state[i] ^= round_keys[offset + i];
    }
}//AddRoundKey
byte_array rijndael::encrypt_block(const byte_array& block) const
{
    if (block.size() != block_size)
    {
        throw std::invalid_argument("invalid block size");
    }
    byte_array state = block;
    add_round_key(state, 0);//хор с ключом раунда 0
    for (size_t round = 1; round < nr; ++round)
    {
        sub_bytes(state);//замена байтов через сбокссс
        shift_rows(state);//цикл сдвиг строк
        mix_columns(state);//умн столбоц на матрицу
        add_round_key(state, round);//хор с ключом тек раунда
    }
    sub_bytes(state);//в последнем только замена и цикл сдвиг
    shift_rows(state);
    add_round_key(state, nr);
    return state;
}
byte_array rijndael::decrypt_block(const byte_array& block) const
{
    if (block.size() != block_size)
    {
        throw std::invalid_argument("invalid block size");
    }
    byte_array state = block;
    add_round_key(state, nr);
    for (size_t round = nr - 1; round > 0; --round)
    {
        inv_shift_rows(state);//обр сдвиг строк
        inv_sub_bytes(state);//обр замена байтов
        add_round_key(state, round);//хор с ключом тек раунда
        inv_mix_columns(state);//обр  умн стобоцов
    }
    inv_shift_rows(state);//обр сдвиг строк,обр замена,хор с ключом тек рауда 0
    inv_sub_bytes(state);
    add_round_key(state, 0);
    return state;
}
const byte_array& rijndael::get_s_box() const
{
    return s_box;
}
const byte_array& rijndael::get_inv_s_box() const
{
    return inv_s_box;
}
size_t rijndael::get_block_size() const
{
    return block_size;
}
size_t rijndael::get_key_size() const
{
    return key_size;
}
size_t rijndael::get_round_count() const
{
    return nr;
}
