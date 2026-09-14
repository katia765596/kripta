#ifndef INTERFACES_H
#define INTERFACES_H
#include <vector>
#include <cstdint>
enum class des_mode { eee3, ede3, eee2, ede2 };
class i_key_schedule {//1.а ген раунд ключей
public:
    virtual std::vector<std::vector<uint8_t>> expand_key(const std::vector<uint8_t>& key) = 0;
    virtual ~i_key_schedule() {}//вирт деструктор для удаления объектов производных классов через указ на интерфейс
};
class i_feistel_round {//1.б
public://метод принимает блок (п.половину) и раунд ключ,возвр преобразованный блок (рез примения функции f)
    virtual std::vector<uint8_t> round_function(const std::vector<uint8_t>& block, const std::vector<uint8_t>& round_key) = 0;
    virtual ~i_feistel_round() {}
};
class i_symmetric_cipher {//1.с 
public:
    virtual void set_key(const std::vector<uint8_t>& key) = 0;
    virtual std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) = 0;
    virtual std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) = 0;
    virtual ~i_symmetric_cipher() {}
};
#endif