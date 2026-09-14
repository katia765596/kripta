#ifndef TRIPLE_DES_H
#define TRIPLE_DES_H
#include "interfaces.h"//1.а,б,с интерф
#include "des.h"//наследник от фейстиляяя
#include <vector>
#include <cstdint>
class triple_des : public i_symmetric_cipher {//наследут интер симм шифр
public:
    triple_des(des_mode mode);//режимки
    void set_key(const std::vector<uint8_t>& key) override;//оверид указывает что метод переопределяет вирт метод базового класса
    std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) override;
    std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) override;
private:
    des_mode mode;
    des des1, des2, des3;
};
#endif
