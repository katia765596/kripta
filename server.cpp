#include <iostream>
#include <thread>
#include <vector>
#include <cstring>
#include <fcntl.h>//для флагов открытия файлов при создании разделяемой памяти
#include <sys/mman.h>//для работы с разделяемой памятью
#include <sys/stat.h>
#include <semaphore.h>//для констант режимов доступа
#include <unistd.h>
#include <signal.h>
#include <cstdlib>
#include "common.h"
static slot* slots = nullptr;//указатель на область разделяемой памяти, хранит массив структур слотс
static sem_t* sem_server[NUM_SLOTS];//массив указ на семафоры, используемые сервером
static sem_t* sem_client[NUM_SLOTS];//массив указ на семафоры для клиентоа
static bool running = true;
void rc4_init(uint8_t* s, const uint8_t* key, size_t key_len) {//иниц s-блока
    for (int i = 0; i < 256; ++i) s[i] = i;
    int j = 0;//перем для перестановок
    for (int i = 0; i < 256; ++i) {
        j = (j + s[i] + key[i % key_len]) & 0xFF;//обн j c учетом текущего байта s-блока и соотв байта ключа
        uint8_t tmp = s[i];//меняются s[i],s[j] в зав от ключа
        s[i] = s[j];
        s[j] = tmp;
    }
}//шифр/дешифр (xor с гаммой)
void rc4_crypt(const uint8_t* key, size_t key_len, const uint8_t* in, uint8_t* out, size_t len) {
    uint8_t s[256];
    rc4_init(s, key, key_len);
    int i = 0, j = 0;
    for (size_t k = 0; k < len; ++k) {
        i = (i + 1) & 0xFF;
        j = (j + s[i]) & 0xFF;
        uint8_t tmp = s[i];
        s[i] = s[j];
        s[j] = tmp;
        out[k] = in[k] ^ s[(s[i] + s[j]) & 0xFF];//в выз буфере берется байт из s-бокса и xor-ится с входным
    }
}
void worker(int id) {//точка входа для каждого потока-обработчика
    while (running) {
        sem_wait(sem_server[id]);//поток блок на семафоре сервера
        if (!running) break;
        slot& s = slots[id];
        if (s.status == 1) {
            rc4_crypt(s.key, s.key_len, s.data, s.result, s.data_len);
            s.result_len = s.data_len;
            s.status = 2;
            sem_post(sem_client[id]);//сигнал клиенту что рез готов
        }
    }
}
void cleanup() {//при заверш сервера для освобожд ресурсов
    running = false;
    for (int i = 0; i < NUM_SLOTS; ++i) {
        if (sem_server[i]) sem_post(sem_server[i]);//будим потоки чтобы они могли выыйти и завершиться
        if (sem_client[i]) sem_post(sem_client[i]);
    }
    if (slots) {
        munmap(slots, sizeof(slot) * NUM_SLOTS);//откл отображение разделяемой памяти
        shm_unlink(SHM_NAME);//удаляет объект раздел памяти из системы
    }
    for (int i = 0; i < NUM_SLOTS; ++i) {
        if (sem_server[i]) {
            sem_close(sem_server[i]);//если сем открыт, то закрываем дескриптом сем и удал из системы
            sem_unlink((std::string(SEM_SERVER_BASE) + std::to_string(i)).c_str());
        }
        if (sem_client[i]) {
            sem_close(sem_client[i]);
            sem_unlink((std::string(SEM_CLIENT_BASE) + std::to_string(i)).c_str());
        }
    }
}
void sig_handler(int) {//принимает номер сигнала и заверашет
    running = false;
}
int main() {
    signal(SIGINT, sig_handler);//сигнал прерывания
    signal(SIGTERM, sig_handler);//сигнал завершения
    int shm_fd = shm_open(SHM_NAME, O_CREAT | O_RDWR, 0666);//создает или открывает объект разделяемой памяти с именем схмнейм, криэей-если не суш, создать, o_rdwr открыть для чтения и запии, 0666 право доступа(чтение,запись)
    if (shm_fd == -1) {
        std::cerr << "shm_open failed\n";
        return 1;
    }
    if (ftruncate(shm_fd, sizeof(slot) * NUM_SLOTS) == -1) {//устанавл размер файла (объекта раздл памяти)
        std::cerr << "ftruncate failed\n";
        return 1;
    }
    slots = (slot*)mmap(0, sizeof(slot) * NUM_SLOTS, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);//отображает разделяемую память в адресное пространство процесса, 0-ядро само выбирает адрес, размер, чтение и запись, изм видны другим процессам, которые тож отображают объект этот, + файл дескр и смещение возвр указ на отобр область
    if (slots == MAP_FAILED) {
        std::cerr << "mmap failed\n";
        return 1;
    }
    close(shm_fd);//дескрипт не нужен, память отображена,закр его
    for (int i = 0; i < NUM_SLOTS; ++i) {
        memset(&slots[i], 0, sizeof(slot));//заполняет нулями всю структуру
        slots[i].status = 0;
        std::string sname = std::string(SEM_SERVER_BASE) + std::to_string(i);//имена сем
        std::string cname = std::string(SEM_CLIENT_BASE) + std::to_string(i);
        sem_unlink(sname.c_str());//удал ранее сущ сем с такими именами
        sem_unlink(cname.c_str());
        sem_server[i] = sem_open(sname.c_str(), O_CREAT | O_EXCL, 0666, 0);
        sem_client[i] = sem_open(cname.c_str(), O_CREAT | O_EXCL, 0666, 0);
        if (sem_server[i] == SEM_FAILED || sem_client[i] == SEM_FAILED) {
            std::cerr << "sem_open failed for slot " << i << "\n";
            cleanup();
            return 1;
        }
    }
    std::vector<std::thread> threads;
    for (int i = 0; i < NUM_SLOTS; ++i) {
        threads.emplace_back(worker, i);// emplace_back конструирует объект theards прямо в векторе, передавая аргументы конструктуру потока
    }
    std::cout << "Server started, waiting for clients...\n";
    while (running) {
        sleep(1);//беск цикл, спит по 1 сек,проверяя статус
    }
    cleanup();
    for (auto& t : threads) t.join();//осн поток ожидает завершения каждого дочернего потока
    std::cout << "Server stopped.\n";
    return 0;
}