#ifndef COMMON_H
#define COMMON_H
#include <cstdint>
#include <cstddef>
#define MAX_DATA_SIZE (1024 * 1024)//макс размер данных
#define MAX_KEY_SIZE 256
#define NUM_SLOTS 100
#define SHM_NAME "/rc4_shm"//имя объекта разделяемой памяти в ос
#define SEM_SERVER_BASE "/rc4_sem_server_"//имя семафора, которое сервер использует для ожидания сигнала от клиента, клиент отправляем сигнал семафору, когда данные записаны
#define SEM_CLIENT_BASE "/rc4_sem_client_"//имя семафора, которое клиент использует для ожидания ответа от сервера, сервер отпр соотв семафор после завершения обработки
struct slot {
    uint8_t key[MAX_KEY_SIZE];
    uint32_t key_len;
    uint8_t data[MAX_DATA_SIZE];
    uint32_t data_len;
    uint8_t result[MAX_DATA_SIZE];
    uint32_t result_len;
    volatile int status; //статут слота - флаг, котор синхронизирует доступ. volatile-указ компилятору, что знач меняется извне (др процессом) и его нельзя кешировать в регисте (читает из опер памяти а не из рег или кэша)
};
#endif