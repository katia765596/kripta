#include "crypto_modes.h"
#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <future>
#include <vector>
namespace
{
    void xor_block(byte_array& first, const byte_array& second)
    {
        for (size_t i = 0; i < first.size(); ++i)
        {
            first[i] ^= second[i];
        }
    }
    byte_array slice_block(
        const byte_array& data,
        size_t position,
        size_t size
    )
    {
        return byte_array(
            data.begin() + position,
            data.begin() + position + size
        );
    }
}
crypto_modes::crypto_modes(
    const block_cipher& algorithm_value,
    crypto_mode mode_value,
    padding_type padding_mode_value,
    const byte_array& initialization_vector_value
)
    : algorithm(algorithm_value),
    mode(mode_value),
    padding_mode(padding_mode_value),
    initialization_vector(initialization_vector_value)
{
    size_t block_size = algorithm.get_block_size();
    if (!initialization_vector.empty() &&
        initialization_vector.size() != block_size)
    {
        throw std::invalid_argument("invalid initialization vector size");
    }
    if (mode != crypto_mode::ecb && initialization_vector.empty())
    {
        initialization_vector.resize(block_size, 0);
    }
}
byte_array crypto_modes::apply_padding(const byte_array& data) const
{
    return padding::add(
        data,
        algorithm.get_block_size(),
        padding_mode
    );
}
byte_array crypto_modes::remove_padding(const byte_array& data) const
{
    return padding::remove(
        data,
        algorithm.get_block_size(),
        padding_mode
    );
}
byte_array crypto_modes::get_iv() const
{
    if (mode == crypto_mode::ecb)
    {
        return byte_array(algorithm.get_block_size(), 0);
    }

    return initialization_vector;
}
byte_array crypto_modes::encrypt(const byte_array& data) const
{
    byte_array padded = apply_padding(data);
    switch (mode)
    {
    case crypto_mode::ecb:
        return encrypt_ecb(padded);
    case crypto_mode::cbc:
        return encrypt_cbc(padded);
    case crypto_mode::pcbc:
        return encrypt_pcbc(padded);
    case crypto_mode::cfb:
        return encrypt_cfb(padded);
    case crypto_mode::ofb:
        return encrypt_ofb(padded);
    case crypto_mode::ctr:
        return crypt_ctr(padded);
    case crypto_mode::random_delta:
        return encrypt_random_delta(padded);
    }

    throw std::invalid_argument("invalid crypto mode");
}
byte_array crypto_modes::decrypt(const byte_array& data) const
{
    byte_array result;
    switch (mode)
    {
    case crypto_mode::ecb:
        result = decrypt_ecb(data);
        break;
    case crypto_mode::cbc:
        result = decrypt_cbc(data);
        break;
    case crypto_mode::pcbc:
        result = decrypt_pcbc(data);
        break;
    case crypto_mode::cfb:
        result = decrypt_cfb(data);
        break;
    case crypto_mode::ofb:
        result = decrypt_ofb(data);
        break;
    case crypto_mode::ctr:
        result = crypt_ctr(data);
        break;
    case crypto_mode::random_delta:
        result = decrypt_random_delta(data);
        break;
    default:
        throw std::invalid_argument("invalid crypto mode");
    }

    return remove_padding(result);
}
byte_array crypto_modes::encrypt_parallel(
    const byte_array& data,
    size_t thread_count
) const
{
    if (thread_count == 0)
    {
        throw std::invalid_argument("invalid thread count");
    }
    if (mode != crypto_mode::ecb && mode != crypto_mode::ctr)
    {
        return encrypt(data);
    }
    byte_array padded = apply_padding(data);
    size_t block_size = algorithm.get_block_size();
    size_t blocks = padded.size() / block_size;
    if (blocks == 0)
    {
        return padded;
    }
    thread_count = std::min(thread_count, blocks);
    byte_array result(padded.size());
    std::vector<std::future<void> > tasks;//дескрипторы потоков
    size_t base = blocks / thread_count;//мин кол-во блоков на поток
    size_t extra = blocks % thread_count;
    size_t first_block = 0;//индекс первого блока для тек потока
    for (size_t thread = 0; thread < thread_count; ++thread)
    {
        size_t count = base + (thread < extra ? 1 : 0);//кол-во блоков которые он должен обработать
        size_t start = first_block;
        first_block += count;
        tasks.push_back(std::async(
            std::launch::async,
            [this, &padded, &result, start, count, block_size]()
            {
                uint64_t counter = static_cast<uint64_t>(start);//для CTR лок счетчик начинающийся с номера первого блока каждый блок испол свой уник счетчик
                for (size_t k = 0; k < count; ++k)
                {
                    size_t position = (start + k) * block_size;//смещение в байтовых массивах
                    byte_array block = slice_block(padded, position, block_size);
                    byte_array encrypted;
                    if (mode == crypto_mode::ecb)
                    {
                        encrypted = algorithm.encrypt_block(block);
                    }
                    else
                    {//если режим CTR, шифрует счетчик получает гамму streamz 
                        byte_array counter_block = make_counter(counter++);
                        byte_array stream = algorithm.encrypt_block(counter_block);
                        encrypted.resize(block_size);
                        for (size_t i = 0; i < block_size; ++i)
                        {//хор блока данных с гаммой
                            encrypted[i] = block[i] ^ stream[i];
                        }
                    }
                    std::copy(
                        encrypted.begin(),
                        encrypted.end(),
                        result.begin() + position
                    );
                }
            }
        ));
    }//закрывает лямбду-функцию и добавляет стд фьюча в вектор tasks
    for (size_t i = 0; i < tasks.size(); ++i)
    {
        tasks[i].get();
    }
    return result;
}
byte_array crypto_modes::decrypt_parallel(
    const byte_array& data,
    size_t thread_count
) const
{
    if (thread_count == 0)
    {
        throw std::invalid_argument("invalid thread count");
    }
    if (mode != crypto_mode::ecb && mode != crypto_mode::ctr)
    {
        return decrypt(data);
    }
    size_t block_size = algorithm.get_block_size();
    if (data.empty() || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid encrypted data size");
    }
    size_t blocks = data.size() / block_size;
    thread_count = std::min(thread_count, blocks);
    byte_array result(data.size());
    std::vector<std::future<void> > tasks;
    size_t base = blocks / thread_count;
    size_t extra = blocks % thread_count;
    size_t first_block = 0;
    for (size_t thread = 0; thread < thread_count; ++thread)
    {
        size_t count = base + (thread < extra ? 1 : 0);
        size_t start = first_block;
        first_block += count;
        tasks.push_back(std::async(
            std::launch::async,
            [this, &data, &result, start, count, block_size]()
            {
                uint64_t counter = static_cast<uint64_t>(start);
                for (size_t k = 0; k < count; ++k)
                {
                    size_t position = (start + k) * block_size;
                    byte_array block = slice_block(data, position, block_size);
                    byte_array decrypted;
                    if (mode == crypto_mode::ecb)
                    {
                        decrypted = algorithm.decrypt_block(block);
                    }
                    else
                    {
                        byte_array counter_block = make_counter(counter++);
                        byte_array stream = algorithm.encrypt_block(counter_block);
                        decrypted.resize(block_size);
                        for (size_t i = 0; i < block_size; ++i)
                        {
                            decrypted[i] = block[i] ^ stream[i];
                        }
                    }
                    std::copy(
                        decrypted.begin(),
                        decrypted.end(),
                        result.begin() + position
                    );
                }
            }
        ));
    }
    for (size_t i = 0; i < tasks.size(); ++i)
    {
        tasks[i].get();
    }
    return remove_padding(result);
}
byte_array crypto_modes::encrypt_ecb(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array block = slice_block(data, i, block_size);
        byte_array encrypted = algorithm.encrypt_block(block);
        std::copy(encrypted.begin(), encrypted.end(), result.begin() + i);
    }
    return result;
}
byte_array crypto_modes::decrypt_ecb(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    if (data.empty() || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid ECB data size");
    }
    byte_array result(data.size());
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array block = slice_block(data, i, block_size);
        byte_array decrypted = algorithm.decrypt_block(block);
        std::copy(decrypted.begin(), decrypted.end(), result.begin() + i);
    }
    return result;
}
byte_array crypto_modes::encrypt_cbc(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array block = slice_block(data, i, block_size);
        xor_block(block, previous);
        byte_array encrypted = algorithm.encrypt_block(block);
        std::copy(encrypted.begin(), encrypted.end(), result.begin() + i);
        previous = encrypted;
    }
    return result;
}
byte_array crypto_modes::decrypt_cbc(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    if (data.empty() || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid CBC data size");
    }
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array block = slice_block(data, i, block_size);
        byte_array decrypted = algorithm.decrypt_block(block);
        xor_block(decrypted, previous);
        std::copy(decrypted.begin(), decrypted.end(), result.begin() + i);
        previous = block;
    }
    return result;
}
byte_array crypto_modes::encrypt_pcbc(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array plain = slice_block(data, i, block_size);
        byte_array input = plain;
        xor_block(input, previous);
        byte_array encrypted = algorithm.encrypt_block(input);
        std::copy(encrypted.begin(), encrypted.end(), result.begin() + i);
        previous = plain;//присваивает открытый текст тек блока
        xor_block(previous, encrypted);//хорит с зашифр тек блоком
    }
    return result;
}
byte_array crypto_modes::decrypt_pcbc(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    if (data.empty() || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid PCBC data size");
    }
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array cipher = slice_block(data, i, block_size);
        byte_array plain = algorithm.decrypt_block(cipher);
        xor_block(plain, previous);
        std::copy(plain.begin(), plain.end(), result.begin() + i);
        previous = plain;
        xor_block(previous, cipher);
    }
    return result;
}
byte_array crypto_modes::encrypt_cfb(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array stream = algorithm.encrypt_block(previous);///штфрует тек сост вектора и получает гамму
        byte_array block = slice_block(data, i, block_size);
        xor_block(block, stream);//хорит с гаммой получая шифротекст
        std::copy(block.begin(), block.end(), result.begin() + i);
        previous = block;//обн вектор этим шифроеткстом
    }
    return result;
}
byte_array crypto_modes::decrypt_cfb(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    if (data.empty() || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid CFB data size");
    }
    byte_array result(data.size());
    byte_array previous = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array cipher = slice_block(data, i, block_size);
        byte_array stream = algorithm.encrypt_block(previous);//ггамма
        byte_array plain = cipher;//копия шифратекста
        xor_block(plain, stream);
        std::copy(plain.begin(), plain.end(), result.begin() + i);
        previous = cipher;
    }
    return result;
}
byte_array crypto_modes::encrypt_ofb(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    byte_array stream = get_iv();
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        stream = algorithm.encrypt_block(stream);//шифрует не вектор а шифр состояние, гамма получается послед шифр предыд знач
        byte_array block = slice_block(data, i, block_size);
        xor_block(block, stream);
        std::copy(block.begin(), block.end(), result.begin() + i);
    }
    return result;
}
byte_array crypto_modes::decrypt_ofb(const byte_array& data) const
{
    return encrypt_ofb(data);
}
byte_array crypto_modes::make_counter(uint64_t value) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array counter = get_iv();

    for (size_t i = 0; i < 8 && i < block_size; ++i)//цикл по последним 8 байтам или меньше если блоксайз меньше 8
    {
        counter[block_size - 1 - i] =
            static_cast<uint8_t>(value >> (i * 8));
    }
    return counter;
}
byte_array crypto_modes::crypt_ctr(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array result(data.size());
    uint64_t counter = 0;
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array counter_block = make_counter(counter++);
        byte_array stream = algorithm.encrypt_block(counter_block);//шифрует блок счётчика получая гамму
        size_t count = std::min(block_size, data.size() - i);//кол-во байт которое нужно обработать в тек блоке (посл может быть не полным)
        for (size_t j = 0; j < count; ++j)
        {
            result[i + j] = data[i + j] ^ stream[j];
        }
    }
    return result;
}
byte_array crypto_modes::encrypt_random_delta(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    byte_array delta(block_size);
    std::random_device device;
    for (size_t i = 0; i < block_size; ++i)
    {
        delta[i] = static_cast<uint8_t>(device());
    }//блок дельта будет добавлен в начало шифротекста
    byte_array result;
    result.reserve(data.size() + block_size);
    result.insert(result.end(), delta.begin(), delta.end());
    for (size_t i = 0; i < data.size(); i += block_size)
    {
        byte_array plain = slice_block(data, i, block_size);
        byte_array input = plain;
        xor_block(input, delta);//хор с копией откр тексьа
        byte_array encrypted = algorithm.encrypt_block(input);
        result.insert(result.end(), encrypted.begin(), encrypted.end());
        delta = plain;
        xor_block(delta, encrypted);
    }
    return result;
}
byte_array crypto_modes::decrypt_random_delta(const byte_array& data) const
{
    size_t block_size = algorithm.get_block_size();
    if (data.size() < block_size || data.size() % block_size != 0)
    {
        throw std::invalid_argument("invalid Random Delta data size");
    }
    byte_array delta(data.begin(), data.begin() + block_size);//извл первый блок как начальную дельту
    byte_array result(data.size() - block_size);
    size_t position = block_size;//начало зашифр блоков
    size_t output = 0;
    while (position < data.size())
    {
        byte_array cipher = slice_block(data, position, block_size);
        byte_array plain = algorithm.decrypt_block(cipher);
        xor_block(plain, delta);
        std::copy(plain.begin(), plain.end(), result.begin() + output);
        delta = plain;
        xor_block(delta, cipher);
        position += block_size;
        output += block_size;
    }
    return result;
}
crypto_mode crypto_modes::get_mode() const
{
    return mode;
}
padding_type crypto_modes::get_padding_type() const
{
    return padding_mode;
}
const byte_array& crypto_modes::get_initialization_vector() const
{
    return initialization_vector;
}
