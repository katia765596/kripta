#ifndef INTERFACES_H
#define INTERFACES_H
#include <vector>
#include <cstdint>
class i_symmetric_cipher {
public:
    virtual void set_key(const std::vector<uint8_t>& key) = 0;
    virtual std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) = 0;
    virtual std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) = 0;
    virtual size_t block_size() const = 0;
    virtual ~i_symmetric_cipher() {}//удал через указ на баз класс
};
class i_cipher_mode {
public:
    virtual void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) = 0;
    virtual std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) = 0;
    virtual bool can_parallel_encrypt() const = 0;
    virtual bool can_parallel_decrypt() const = 0;
    virtual bool needs_padding() const = 0; 
    virtual void reset() = 0;//сброс внутр сост режима,чтобы переиспол объект для новой операции
    virtual ~i_cipher_mode() {}
};
class i_padding {
public:
    virtual std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data, size_t block_size) = 0;
    virtual std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data, size_t block_size) = 0;
    virtual ~i_padding() {}
};
#endif