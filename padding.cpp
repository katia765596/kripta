#include "padding.h"
#include <random>
#include <stdexcept>
byte_array padding::add(
    const byte_array& data,
    size_t block_size,
    padding_type type
)
{
    if (block_size == 0 || block_size > 255)
    {//посл байт хранит кол-во добавленных байт от 0 до 255
        throw std::invalid_argument("invalid block size");
    }
    switch (type)
    {
    case padding_type::zeros:
        return add_zeros(data, block_size);
    case padding_type::pkcs7:
        return add_pkcs7(data, block_size);
    case padding_type::iso_10126:
        return add_iso_10126(data, block_size);
    case padding_type::ansi_x923:
        return add_ansi_x923(data, block_size);
    }
    throw std::invalid_argument("invalid padding type");
}
byte_array padding::remove(
    const byte_array& data,
    size_t block_size,
    padding_type type
)
{
    if (block_size == 0 || block_size > 255)
    {
        throw std::invalid_argument("invalid block size");
    }
    switch (type)
    {
    case padding_type::zeros:
        return remove_zeros(data);
    case padding_type::pkcs7:
        return remove_pkcs7(data, block_size);
    case padding_type::iso_10126:
        return remove_iso_10126(data, block_size);
    case padding_type::ansi_x923:
        return remove_ansi_x923(data, block_size);
    }

    throw std::invalid_argument("invalid padding type");
}
byte_array padding::add_zeros(
    const byte_array& data,
    size_t block_size
)
{
    byte_array result = data;//создает копию входных данных в резалт
    size_t remainder = result.size() % block_size;
    if (remainder != 0)
    {
        result.insert(result.end(), block_size - remainder, 0);//инсерт вставляет кол-во добавляемых байт блок_сайз-остаток значение 0 в коонец выектора
    }

    return result;
}
byte_array padding::add_pkcs7(//добавляет байты, каждый из которых равен кол-ву добавл байт
    const byte_array& data,
    size_t block_size
)
{
    byte_array result = data;
    size_t count = block_size - result.size() % block_size;
    if (count == 0)
    {
        count = block_size;
    }
    result.insert(result.end(), count, static_cast<uint8_t>(count));//вставляет в конец каунт байт, каждый каунт
    return result;
}
byte_array padding::add_iso_10126(//последний содержит кол-во добавленных байт, остальные добавленные байты случайны
    const byte_array& data,
    size_t block_size
)
{
    byte_array result = data;
    size_t count = block_size - result.size() % block_size;
    if (count == 0)
    {
        count = block_size;
    }
    std::random_device device;
    for (size_t i = 1; i < count; ++i)
    {
        result.push_back(static_cast<uint8_t>(device()));
    }
    result.push_back(static_cast<uint8_t>(count));
    return result;
}
byte_array padding::add_ansi_x923(
    const byte_array& data,
    size_t block_size
)
{
    byte_array result = data;
    size_t count = block_size - result.size() % block_size;
    if (count == 0)
    {
        count = block_size;
    }
    if (count > 1)
    {
        result.insert(result.end(), count - 1, 0);//если добавл более одного байта вставляем count-1 байт со знач 0, все добавляем байты, кроме последнего 0
    }
    result.push_back(static_cast<uint8_t>(count));
    return result;
}
byte_array padding::remove_zeros(const byte_array& data)
{
    byte_array result = data;

    while (!result.empty() && result.back() == 0)
    {
        result.pop_back();
    }

    return result;
}
byte_array padding::remove_pkcs7(
    const byte_array& data,
    size_t block_size
)
{
    if (data.empty())
    {
        throw std::invalid_argument("empty padded data");
    }
    uint8_t count = data.back();//кол-во добавленных байт
    if (count == 0 || count > block_size || count > data.size())
    {
        throw std::invalid_argument("invalid pkcs7 padding");
    }
    for (size_t i = 0; i < count; ++i)
    {
        if (data[data.size() - 1 - i] != count)//проверяем все удал байты дейст равны каунт
        {
            throw std::invalid_argument("invalid pkcs7 padding");
        }
    }
    return byte_array(data.begin(), data.end() - count);
}
byte_array padding::remove_iso_10126(
    const byte_array& data,
    size_t block_size
)
{
    if (data.empty())
    {
        throw std::invalid_argument("empty padded data");
    }
    uint8_t count = data.back();
    if (count == 0 || count > block_size || count > data.size())
    {
        throw std::invalid_argument("invalid iso 10126 padding");
    }
    return byte_array(data.begin(), data.end() - count);
}
byte_array padding::remove_ansi_x923(
    const byte_array& data,
    size_t block_size
)
{
    if (data.empty())
    {
        throw std::invalid_argument("empty padded data");
    }
    uint8_t count = data.back();
    if (count == 0 || count > block_size || count > data.size())
    {
        throw std::invalid_argument("invalid ansi x923 padding");
    }
    for (size_t i = 1; i < count; ++i)
    {
        if (data[data.size() - 1 - i] != 0)
        {
            throw std::invalid_argument("invalid ansi x923 padding");
        }
    }
    return byte_array(data.begin(), data.end() - count);
}
