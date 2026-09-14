#ifndef MODES_H
#define MODES_H
#include "interfaces.h"
#include <vector>
#include <cstdint>
class ecb_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> iv;
};
class cbc_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> prev_block;
    bool initialized;
};
class pcbc_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> prev_cipher;
    std::vector<uint8_t> prev_plain;
    bool initialized;
};
class cfb_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> shift_reg;
};
class ofb_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> state;
};
class ctr_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> counter;
    uint64_t counter_value;
};
class random_delta_mode : public i_cipher_mode {
public:
    void init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) override;
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) override;
    bool can_parallel_encrypt() const override;
    bool can_parallel_decrypt() const override;
    bool needs_padding() const override;
    void reset() override;
private:
    std::vector<uint8_t> initial;
    std::vector<uint8_t> delta;
    uint64_t current_counter;
    bool initialized;
};
#endif