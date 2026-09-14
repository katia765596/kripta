#ifndef BLOCK_CIPHER_H
#define BLOCK_CIPHER_H
#include "byte_array.h"
#include <cstddef>
class block_cipher//интерфейс любого блочного шифра
{
public:
    virtual ~block_cipher() = default;
    virtual byte_array encrypt_block(const byte_array& block) const = 0;
    virtual byte_array decrypt_block(const byte_array& block) const = 0;
    virtual size_t get_block_size() const = 0;
};
#endif
