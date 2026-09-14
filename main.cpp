#include "cipher_context.h"
#include <iostream>
#include <cassert>
#include <future>
#include <fstream>
#include <cstring>
#include <iomanip>
class simple_cipher : public i_symmetric_cipher {
public:
    void set_key(const std::vector<uint8_t>& key) override {
        this->key = key;
    }
    std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) override {
        return xor_block(block);
    }
    std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) override {
        return xor_block(block);
    }
    size_t block_size() const override {
        return 8;
    }
private:
    std::vector<uint8_t> key;
    std::vector<uint8_t> xor_block(const std::vector<uint8_t>& block) {
        size_t bs = block_size();//прив метод который XORит каждый байт блока с соотв байтом ключа (по модулю длины ключа)
        std::vector<uint8_t> result(bs);
        for (size_t i = 0; i < bs; ++i) {
            result[i] = block[i] ^ (key.empty() ? 0x00 : key[i % key.size()]);//делим ключ если он короче блока
        }
        return result;
    }
};
static void test_encrypt_decrypt(cipher_context& ctx, const std::vector<uint8_t>& plain) {
    std::vector<uint8_t> ciphertext, decrypted;
    ctx.encrypt(plain, ciphertext, 1);
    ctx.decrypt(ciphertext, decrypted, 1);
    if (plain != decrypted) {
        std::cerr << "error, plain.size=" << plain.size() << ", decrypted.size=" << decrypted.size() << std::endl;
        size_t n = std::min(plain.size(), decrypted.size());
        for (size_t i = 0; i < n; ++i) {
            if (plain[i] != decrypted[i]) {
                std::cerr << "byte " << i << ": plain=0x" << std::hex << (int)plain[i] << ", dec=0x" << (int)decrypted[i] << std::dec << std::endl;
            }
        }
        if (plain.size() != decrypted.size()) {
            std::cerr << "sizes differ" << std::endl;
        }
        assert(plain == decrypted);
    }
}
int main() {
    simple_cipher cipher;
    std::vector<uint8_t> key(8, 0x55);
    cipher.set_key(key);
    std::vector<uint8_t> iv(8, 0x00);
    std::vector<cipher_mode> modes = {
        cipher_mode::ecb,
        cipher_mode::cbc,
        cipher_mode::pcbc,
        cipher_mode::cfb,
        cipher_mode::ofb,
        cipher_mode::ctr,
        cipher_mode::random_delta
    };
    for (auto m : modes) {
        std::vector<uint8_t> plain;
        if (m == cipher_mode::ecb || m == cipher_mode::cbc || m == cipher_mode::pcbc) {
            plain = { 0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A };
        }
        else {
            plain = { 0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08 };
        }
        if (m == cipher_mode::random_delta) {
            cipher_context ctx(&cipher, m, padding_mode::pkcs7, iv,
                { 0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08 });
            test_encrypt_decrypt(ctx, plain);
        }
        else {
            cipher_context ctx(&cipher, m, padding_mode::pkcs7, iv);
            test_encrypt_decrypt(ctx, plain);
        }
        std::cout << "Mode " << static_cast<int>(m) << " OK\n";
    }
    std::vector<padding_mode> pads = {
        padding_mode::zeros,
        padding_mode::ansi_x923,
        padding_mode::pkcs7,
        padding_mode::iso10126
    };
    for (auto p : pads) {
        cipher_context ctx(&cipher, cipher_mode::cbc, p, iv);
        std::vector<uint8_t> plain = { 0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09 };
        test_encrypt_decrypt(ctx, plain);
        std::cout << "Padding " << static_cast<int>(p) << " OK\n";
    }
    {
        std::ofstream("test_in.txt") << "hello world";
        cipher_context ctx(&cipher, cipher_mode::cbc, padding_mode::pkcs7, iv);
        auto fut = ctx.encrypt_async("test_in.txt", "test_enc.bin", 1);
        fut.wait();//ждем завршения
        ctx.decrypt_async("test_enc.bin", "test_dec.txt", 1).wait();
        std::ifstream fin("test_dec.txt", std::ios::binary);
        fin.seekg(0, std::ios::end);
        size_t size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        std::string content(size, ' ');
        fin.read(&content[0], size);
        fin.close();
        std::cout << "test_dec.txt size: " << size << ", content: '" << content << "'" << std::endl;
        if (content != "hello world") {
            std::cerr << "error" << std::endl;
            std::cerr << "got: " << content << std::endl;
        }
        assert(content == "hello world");
        std::cout << "File async test OK\n";
    }
    std::remove("test_in.txt");
    std::remove("test_enc.bin");
    std::remove("test_dec.txt");
    std::cout << "All tests passed\n";
    return 0;
}