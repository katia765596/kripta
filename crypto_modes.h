#ifndef CRYPTO_MODES_H
#define CRYPTO_MODES_H
#include "byte_array.h"
#include "padding.h"
#include "block_cipher.h"
#include <cstddef>
enum class crypto_mode
{
    ecb,cbc,pcbc,cfb,ofb,ctr,random_delta
};
class crypto_modes
{
public:
    crypto_modes(
        const block_cipher& algorithm,//ссылка на объект блочного шифра
        crypto_mode mode,//режим шифрования
        padding_type padding_mode,
        const byte_array& initialization_vector = byte_array()
    );
    byte_array encrypt(const byte_array& data) const;
    byte_array decrypt(const byte_array& data) const;
    byte_array encrypt_parallel(const byte_array& data, size_t thread_count) const;
    byte_array decrypt_parallel(const byte_array& data, size_t thread_count) const;
    crypto_mode get_mode() const;
    padding_type get_padding_type() const;
    const byte_array& get_initialization_vector() const;
private:
    const block_cipher& algorithm;
    crypto_mode mode;
    padding_type padding_mode;
    byte_array initialization_vector;
    byte_array encrypt_ecb(const byte_array& data) const;
    byte_array decrypt_ecb(const byte_array& data) const;
    byte_array encrypt_cbc(const byte_array& data) const;
    byte_array decrypt_cbc(const byte_array& data) const;
    byte_array encrypt_pcbc(const byte_array& data) const;
    byte_array decrypt_pcbc(const byte_array& data) const;
    byte_array encrypt_cfb(const byte_array& data) const;
    byte_array decrypt_cfb(const byte_array& data) const;
    byte_array encrypt_ofb(const byte_array& data) const;
    byte_array decrypt_ofb(const byte_array& data) const;
    byte_array crypt_ctr(const byte_array& data) const;
    byte_array encrypt_random_delta(const byte_array& data) const;
    byte_array decrypt_random_delta(const byte_array& data) const;
    byte_array apply_padding(const byte_array& data) const;
    byte_array remove_padding(const byte_array& data) const;
    byte_array get_iv() const;
    byte_array make_counter(uint64_t value) const;//блок четчика для ctr
};
#endif
