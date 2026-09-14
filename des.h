#ifndef DES_H
#define DES_H
#include "feistel_network.h"//объявл унив класса сети фейстеля от него наследуем
#include "interfaces.h"
#include <memory>
#include <vector>
#include <cstdint>
class des : public feistel_wrapper {//наследует от обертки с шабл методом которая сама реализует интерфейс симм шифр, так дес становится готовым симм шифром, который испол сеть фейстеля с конкретными параметрами дес
public:
    des();//вызывает констр баз класса  передавая ему реализации интер расширения ключа и раунд ф
protected:
    std::vector<uint8_t> pre_processing(const std::vector<uint8_t>& block) override;//переопр вирт методы перед шифр сетью ф (нач перест ip)
    std::vector<uint8_t> post_processing(const std::vector<uint8_t>& block) override;//после сети фейстеля (конечная перестановка fp)
private:
    class des_key_schedule : public i_key_schedule {//наследует кей шедул от класса генер раунд ключей и реалз метод
    public:
        std::vector<std::vector<uint8_t>> expand_key(const std::vector<uint8_t>& key) override;//из 8 байт 16 раунд ключей по 6 байт
    private:
        uint64_t bytes_to_uint64(const std::vector<uint8_t>& bytes);
        std::vector<uint8_t> uint64_to_bytes(uint64_t val, int len);//48 бит обр преобр
        uint64_t permute_64(uint64_t val, const int* table, int bits);//общ ф перестановки битов по табл
    };//методы приватные нужны внутри класса
    class des_feistel_round : public i_feistel_round {
    public://из 4 байт блока (п.половина) и 6-байтового раунд ключа делает 4байтовый блок правый (рез ф)
        std::vector<uint8_t> round_function(const std::vector<uint8_t>& block, const std::vector<uint8_t>& round_key) override;
    private:
        uint32_t bytes_to_uint32(const std::vector<uint8_t>& bytes);//для половин в 32 бита
        uint64_t bytes_to_uint48(const std::vector<uint8_t>& bytes);//преобр раунд ключа в 48 бит
        std::vector<uint8_t> uint32_to_bytes(uint32_t val);
        uint64_t expand_32_to_48(uint32_t val);//расширение п половины с 32 до 48(е-расширение)
        uint32_t permute_32(uint32_t val, const int* table, int bits);//перестановка 32 бит знач по табл
    };
    std::unique_ptr<i_key_schedule> ks_holder;//умн указ на объекты влож классов (хранят реал интер И ПЕРЕДАЮТ в баз класс который жет указ)
    std::unique_ptr<i_feistel_round> fr_holder;
    uint64_t bytes_to_uint64(const std::vector<uint8_t>& bytes);
    std::vector<uint8_t> uint64_to_bytes(uint64_t val);
    uint64_t permute_64(uint64_t val, const int* table, int bits);
};
#endif