#ifndef LAN_KEY_EXCHANGE_H
#define LAN_KEY_EXCHANGE_H
#include "byte_array.h"
#include "diffie_hellman.h"//реализован протокол диффи-хеллмана, испол для генерации общих секретов
#include <cstdint>
#include <string>
class lan_key_exchange//для организации обмена ключами по сети (лок сеть)
{
public:
    static byte_array run_server(const diffie_hellman& participant, uint16_t port, size_t key_size);//(запускает серверную часть обмена ключами,принмает партисипент(параметры DH, публичный ключ), номер порта, желаемый размер ключа, возвр общий секрет
    static byte_array run_client(const diffie_hellman& participant, const std::string& address, uint16_t port, size_t key_size);//принимает партисипент, строку с адресом сервера,порт, размер ключа, устанавливает соед, обмен ключами, возв общий секрет
};
#endif
