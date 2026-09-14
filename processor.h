#ifndef PROCESSOR_H
#define PROCESSOR_H
#include "interfaces.h"
#include <vector>
#include <cstdint>
#include <string>
class cipher_processor {//обертка над любым симм шифром (класс использует указатель на симметрик шифр), через констр передается реал симм алг шифр
public://класс обрабатывает данные с помощью переданного алг шифр
    cipher_processor(i_symmetric_cipher* cipher);//указ на объект реалз интерфейс классу передаетс конкретный шифр
    std::vector<uint8_t> process_data(const std::vector<uint8_t>& data, bool encrypt);
    void process_file(const std::string& in_path, const std::string& out_path, bool encrypt);
private:
    i_symmetric_cipher* algo;//указ на объект алг шифр
    size_t block_size = 8;
    std::vector<uint8_t> add_padding(const std::vector<uint8_t>& data);
    std::vector<uint8_t> remove_padding(const std::vector<uint8_t>& data);
};
bool compare_files(const std::string& f1, const std::string& f2);
#endif