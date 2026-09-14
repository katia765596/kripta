#ifndef PADDING_H
#define PADDING_H
#include "byte_array.h"
#include <cstddef>
enum class padding_type
{
    zeros,
    pkcs7,
    iso_10126,
    ansi_x923
};
class padding
{
public:
    static byte_array add(const byte_array& data, size_t block_size, padding_type type);
    static byte_array remove(const byte_array& data, size_t block_size, padding_type type);
private:
    static byte_array add_zeros(const byte_array& data, size_t block_size);
    static byte_array add_pkcs7(const byte_array& data, size_t block_size);
    static byte_array add_iso_10126(const byte_array& data, size_t block_size);
    static byte_array add_ansi_x923(const byte_array& data, size_t block_size);
    static byte_array remove_zeros(const byte_array& data);
    static byte_array remove_pkcs7(const byte_array& data, size_t block_size);
    static byte_array remove_iso_10126(const byte_array& data, size_t block_size);
    static byte_array remove_ansi_x923(const byte_array& data, size_t block_size);
};
#endif
