#include "paddings.h"
#include <stdexcept>
#include <cstdlib>
#include <ctime>
std::vector<uint8_t> zeros_padding::add_padding(const std::vector<uint8_t>& data, size_t block_size) {
    size_t pad_len = block_size - (data.size() % block_size);
    if (pad_len == 0) pad_len = block_size;
    std::vector<uint8_t> padded = data;
    padded.resize(data.size() + pad_len, 0x00);//знач которым заполняются новые байты (нулями)
    return padded;
}
std::vector<uint8_t> zeros_padding::remove_padding(const std::vector<uint8_t>& data, size_t block_size) {
    if (data.empty()) return data;
    size_t i = data.size();//запоминаем тек размер
    while (i > 0 && data[i - 1] == 0x00) --i;//пока i>0 и пред равен нулю умен i идем с конца и срезаем нулевые байты
    return std::vector<uint8_t>(data.begin(), data.begin() + i);//до поз i не включая
}
std::vector<uint8_t> ansi_x923_padding::add_padding(const std::vector<uint8_t>& data, size_t block_size) {
    size_t pad_len = block_size - (data.size() % block_size);
    if (pad_len == 0) pad_len = block_size;
    std::vector<uint8_t> padded = data;
    padded.resize(data.size() + pad_len, 0x00);
    padded.back() = static_cast<uint8_t>(pad_len);//послед байт устанавливаем равным пад_лен- кол-ву добавленных байт
    return padded;
}
std::vector<uint8_t> ansi_x923_padding::remove_padding(const std::vector<uint8_t>& data, size_t block_size) {
    if (data.empty()) throw std::runtime_error("empty data");
    uint8_t pad_len = data.back();//посл байт - кол-во добавл байт
    if (pad_len == 0 || pad_len > block_size) throw std::runtime_error("invalid padding");
    for (size_t i = data.size() - pad_len; i < data.size() - 1; ++i) {//проходим по байтам, которые должно быть нулевыми(кроме последнего)
        if (data[i] != 0x00) throw std::runtime_error("invalid padding bytes");
    }
    return std::vector<uint8_t>(data.begin(), data.end() - pad_len);//возвр вектор без пад_лен байт
}
std::vector<uint8_t> pkcs7_padding::add_padding(const std::vector<uint8_t>& data, size_t block_size) {
    size_t pad_len = block_size - (data.size() % block_size);
    if (pad_len == 0) pad_len = block_size;
    std::vector<uint8_t> padded = data;
    padded.resize(data.size() + pad_len, static_cast<uint8_t>(pad_len));//здесь при рисайз мы указываем пад_лен для новых байтов все добавленные = пад_лен
    return padded;
}
std::vector<uint8_t> pkcs7_padding::remove_padding(const std::vector<uint8_t>& data, size_t block_size) {
    if (data.empty()) throw std::runtime_error("empty data");
    uint8_t pad_len = data.back();
    if (pad_len == 0 || pad_len > block_size) throw std::runtime_error("invalid padding");
    for (size_t i = data.size() - pad_len; i < data.size(); ++i) {
        if (data[i] != pad_len) throw std::runtime_error("invalid padding bytes");//все байты в диапазоне паддинга пад лен
    }
    return std::vector<uint8_t>(data.begin(), data.end() - pad_len);
}
uint8_t iso10126_padding::random_byte() {//прив метод для генерации случ байта
    static bool seeded = false;//стат переменная запоминает был ли уже иниц генератор сохр сост между вызовами
    if (!seeded) { std::srand(std::time(nullptr)); seeded = true; }//если нет то вызываем с тек временом
    return static_cast<uint8_t>(std::rand() % 256);
}
std::vector<uint8_t> iso10126_padding::add_padding(const std::vector<uint8_t>& data, size_t block_size) {
    size_t pad_len = block_size - (data.size() % block_size);
    if (pad_len == 0) pad_len = block_size;
    std::vector<uint8_t> padded = data;
    padded.resize(data.size() + pad_len);
    for (size_t i = data.size(); i < padded.size() - 1; ++i) {
        padded[i] = random_byte();//расширили вектор проходим по байтам кроме последнего и заполняем пад_лен случ числами
    }
    padded.back() = static_cast<uint8_t>(pad_len);//послд байт = пад лен
    return padded;
}
std::vector<uint8_t> iso10126_padding::remove_padding(const std::vector<uint8_t>& data, size_t block_size) {
    if (data.empty()) throw std::runtime_error("empty data");
    uint8_t pad_len = data.back();
    if (pad_len == 0 || pad_len > block_size) throw std::runtime_error("invalid padding");
    return std::vector<uint8_t>(data.begin(), data.end() - pad_len);
}