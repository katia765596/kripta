#include "processor.h"//пункт 5
#include <fstream>
#include <stdexcept>
#include <cstring>
cipher_processor::cipher_processor(i_symmetric_cipher* cipher) : algo(cipher) {}//иниц поле алго указателем на объект алг(интерфейс_
std::vector<uint8_t> cipher_processor::process_data(const std::vector<uint8_t>& data, bool encrypt) {
    if (!algo) throw std::runtime_error("no cipher set");
    if (data.empty()) return std::vector<uint8_t>();
    std::vector<uint8_t> result;
    if (encrypt) {
        std::vector<uint8_t> padded = add_padding(data);
        result.reserve(padded.size());//резерв памяти для рез с паддингом
        for (size_t i = 0; i < padded.size(); i += block_size) {//поблочно 8 байт
            std::vector<uint8_t> block(padded.begin() + i, padded.begin() + i + block_size);
            std::vector<uint8_t> enc = algo->encrypt_block(block);//вызывает метод энкрипт переданного алг, передавая тек блок
            result.insert(result.end(), enc.begin(), enc.end());//вставляем в конец вектора
        }
    }
    else {
        if (data.size() % block_size != 0) throw std::invalid_argument("ciphertext size not multiple of block");
        std::vector<uint8_t> decrypted;
        decrypted.reserve(data.size());
        for (size_t i = 0; i < data.size(); i += block_size) {
            std::vector<uint8_t> block(data.begin() + i, data.begin() + i + block_size);
            std::vector<uint8_t> dec = algo->decrypt_block(block);
            decrypted.insert(decrypted.end(), dec.begin(), dec.end());
        }
        result = remove_padding(decrypted);//потом удал паддинг
    }
    return result;
}
void cipher_processor::process_file(const std::string& in_path, const std::string& out_path, bool encrypt) {
    std::ifstream fin(in_path, std::ios::binary);
    if (!fin) throw std::runtime_error("cannot open input file");
    std::ofstream fout(out_path, std::ios::binary);
    if (!fout) throw std::runtime_error("cannot open output file");
    if (encrypt) {
        std::vector<uint8_t> buffer(block_size * 1024);//читаем файл кусками
        std::vector<uint8_t> last_block;//неполный блок меньший 8 байт добавим паддинг
        while (fin) {
            fin.read(reinterpret_cast<char*>(buffer.data()), buffer.size());//преобр в чар
            size_t read_bytes = fin.gcount();
            if (read_bytes == 0) break;
            size_t blocks_to_process = read_bytes / block_size;//rол-во полных блоков
            size_t leftover = read_bytes % block_size;//неполный блок
            if (leftover > 0) {
                last_block.insert(last_block.end(), buffer.begin() + blocks_to_process * block_size,
                    buffer.begin() + blocks_to_process * block_size + leftover);// добавляет в ласт блок остаток
            }
            for (size_t i = 0; i < blocks_to_process * block_size; i += block_size) {
                std::vector<uint8_t> block(buffer.begin() + i, buffer.begin() + i + block_size);
                std::vector<uint8_t> enc = algo->encrypt_block(block);//проходит по всем блокам, шифрует и записывает в вых ф)
                fout.write(reinterpret_cast<const char*>(enc.data()), enc.size());
            }
        }
        if (!last_block.empty() || fin.eof()) {//падд даже если кратно 8
            std::vector<uint8_t> padded = add_padding(last_block);
            for (size_t i = 0; i < padded.size(); i += block_size) {
                std::vector<uint8_t> block(padded.begin() + i, padded.begin() + i + block_size);
                std::vector<uint8_t> enc = algo->encrypt_block(block);
                fout.write(reinterpret_cast<const char*>(enc.data()), enc.size());
            }
        }
    }
    else { 
        std::vector<uint8_t> buffer(block_size * 1024);
        std::vector<uint8_t> prev_block;
        bool last_block_processed = false;//флаг будет ли обработан посл блок
        while (fin) {
            fin.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
            size_t read_bytes = fin.gcount();
            if (read_bytes == 0) break;
            size_t total_blocks = read_bytes / block_size;
            size_t bytes_to_process = read_bytes - (read_bytes % block_size);//полные блоки минус остаток
            if (fin.eof()) {
                if (read_bytes < block_size) throw std::runtime_error("invalid ciphertext");
                for (size_t i = 0; i < bytes_to_process - block_size; i += block_size) {
                    std::vector<uint8_t> block(buffer.begin() + i, buffer.begin() + i + block_size);
                    std::vector<uint8_t> dec = algo->decrypt_block(block);
                    fout.write(reinterpret_cast<const char*>(dec.data()), dec.size());
                }
                std::vector<uint8_t> last_block_vec(buffer.begin() + bytes_to_process - block_size,
                    buffer.begin() + bytes_to_process);//извлекает последние 8 байт из буфера
                std::vector<uint8_t> dec_last = algo->decrypt_block(last_block_vec);
                std::vector<uint8_t> unpadded = remove_padding(dec_last);
                fout.write(reinterpret_cast<const char*>(unpadded.data()), unpadded.size());//запись данных без паддинга
                last_block_processed = true;
            }
            else {//если не последний дешефруем все блоки из буфера и сразу запись
                for (size_t i = 0; i < bytes_to_process; i += block_size) {
                    std::vector<uint8_t> block(buffer.begin() + i, buffer.begin() + i + block_size);
                    std::vector<uint8_t> dec = algo->decrypt_block(block);
                    fout.write(reinterpret_cast<const char*>(dec.data()), dec.size());
                }
            }
        }
        if (!last_block_processed && fin.eof()) {
            throw std::runtime_error("no data processed");
        }
    }
    fin.close();
    fout.close();
}
std::vector<uint8_t> cipher_processor::add_padding(const std::vector<uint8_t>& data) {
    size_t pad_len = block_size - (data.size() % block_size);//скок байт добавить
    if (pad_len == 0) pad_len = block_size;
    std::vector<uint8_t> padded = data;//коп исх данные
    padded.resize(data.size() + pad_len, static_cast<uint8_t>(pad_len));//изм размер вектора
    return padded;
}
std::vector<uint8_t> cipher_processor::remove_padding(const std::vector<uint8_t>& data) {
    if (data.empty()) return data;
    uint8_t pad_val = data.back();//получаем посл байт - знач паддинга
    if (pad_val == 0 || pad_val > block_size) throw std::runtime_error("invalid padding");
    size_t pad_len = pad_val;//преобр в сайз т
    if (pad_len > data.size()) throw std::runtime_error("invalid padding length");
    for (size_t i = data.size() - pad_len; i < data.size(); ++i) {
        if (data[i] != pad_val) throw std::runtime_error("invalid padding bytes");
    }//проверка все байты паддинга: каждый должен быть равен пад_вал
    return std::vector<uint8_t>(data.begin(), data.end() - pad_len);
}//возвр вектор без паддинга
bool compare_files(const std::string& f1, const std::string& f2) {
    std::ifstream in1(f1, std::ios::binary);
    std::ifstream in2(f2, std::ios::binary);
    if (!in1 || !in2) return false;
    const size_t buf_size = 8192;//размер буфера для чтения по частям
    char buf1[buf_size], buf2[buf_size];
    while (in1.good() && in2.good()) {
        in1.read(buf1, buf_size);
        in2.read(buf2, buf_size);//читаем по буф сайз байт
        if (in1.gcount() != in2.gcount()) return false;
        if (std::memcmp(buf1, buf2, in1.gcount()) != 0) return false;
    }
    return in1.eof() && in2.eof();
}