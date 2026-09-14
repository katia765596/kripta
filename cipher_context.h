#ifndef CIPHER_CONTEXT_H
#define CIPHER_CONTEXT_H
#include "interfaces.h"
#include "enums.h"
#include <vector>
#include <cstdint>
#include <string>
#include <future>
#include <memory>
class cipher_context {
public:
    cipher_context(i_symmetric_cipher* algo,
        cipher_mode mode,
        padding_mode pad,
        const std::vector<uint8_t>& iv = {},
        const std::initializer_list<uint8_t>& extra = {});
    ~cipher_context();//только для очистки умных указ не для автом упр памятью для шифер
    void encrypt(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, int num_threads = 1);
    void decrypt(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, int num_threads = 1);
    std::future<void> encrypt_async(const std::string& input_path, const std::string& output_path, int num_threads = 1);//возвр фюча чтоб запустить опер в отдельном потоке и дождаться завершения
    std::future<void> decrypt_async(const std::string& input_path, const std::string& output_path, int num_threads = 1);//внутри методы испол асунк для запуска прив метода проц_файл_сунк в фоновом потоке
private:
    i_symmetric_cipher* cipher;
    std::unique_ptr<i_cipher_mode> mode;
    std::unique_ptr<i_padding> padding;
    std::vector<uint8_t> iv;
    std::vector<uint8_t> extra_params;
    size_t block_size;
    void process_blocks(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, bool encrypt, int num_threads);
    void process_file_sync(const std::string& input_path, const std::string& output_path, bool encrypt, int num_threads);
};
#endif