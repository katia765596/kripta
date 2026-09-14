#ifndef FEISTEL_NETWORK_H//обертка над сетью фейстеля пункт 2
#define FEISTEL_NETWORK_H
#include "interfaces.h"
#include <vector>
#include <cstdint>
class feistel_network {
public:
    feistel_network(i_key_schedule* ks, i_feistel_round* fr, int rounds);
    std::vector<uint8_t> process_block(const std::vector<uint8_t>& block, bool encrypt, const std::vector<uint8_t>& key);//генер раунд ключей внутри,+ ф к п половине
private:
    i_key_schedule* key_sched;//указ на генератор ключей
    i_feistel_round* round_func;//указ на раунд ф
    int num_rounds;//число раундов
};
class feistel_wrapper : public i_symmetric_cipher {//класс-наследник реализует интерфейс шифра 1.с
public:
    feistel_wrapper(i_key_schedule* ks, i_feistel_round* fr, int rounds);
    void set_key(const std::vector<uint8_t>& key) override; //оверид указывает что метод определяет вирт метод баз класса
    std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) override;//реал шабл метод
    std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) override;
protected://члены класса доступны внутри самого класса и в производных
    virtual std::vector<uint8_t> pre_processing(const std::vector<uint8_t>& block);//вирт чтобы наследники могли переопд.-паттерн шабл метод
    virtual std::vector<uint8_t> post_processing(const std::vector<uint8_t>& block);
    feistel_network network;
    std::vector<uint8_t> current_key;
};
#endif