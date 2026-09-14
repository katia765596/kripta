#ifndef PADDINGS_H
#define PADDINGS_H
#include "interfaces.h"
#include <vector>
#include <cstdint>
class zeros_padding : public i_padding {
public:
    std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data, size_t block_size) override;
    std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data, size_t block_size) override;
};
class ansi_x923_padding : public i_padding {
public:
    std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data, size_t block_size) override;
    std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data, size_t block_size) override;
};
class pkcs7_padding : public i_padding {
public:
    std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data, size_t block_size) override;
    std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data, size_t block_size) override;
};
class iso10126_padding : public i_padding {
public:
    std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data, size_t block_size) override;
    std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data, size_t block_size) override;
private:
    uint8_t random_byte();
};
#endif