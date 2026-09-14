#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <semaphore.h>
#include <unistd.h>
#include <cstdlib>
#include "common.h"
static slot* slots = nullptr;
static sem_t* sem_server = nullptr;
static sem_t* sem_client = nullptr;
bool compare_files(const std::string& f1, const std::string& f2) {
    std::ifstream in1(f1, std::ios::binary | std::ios::ate);//создаются два объекта входных файловых потоков
    std::ifstream in2(f2, std::ios::binary | std::ios::ate);//открыв в бин режиме,устанавливаем позицию чтения в конец файла (чтобы потом получить размер)
    if (!in1 || !in2) return false;
    if (in1.tellg() != in2.tellg()) return false;//телг возвр тек поз (в байтах) отерыли с конца
    in1.seekg(0);//перемещаем поз чтения в начало
    in2.seekg(0);
    char c1, c2;
    while (in1.get(c1) && in2.get(c2)) {//читает один байт из 1 ф в c1 и возвр ссылку на поток, который в усл while преобразцется в bool
        if (c1 != c2) return false;
    }
    return true;
}
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage:\n";
        std::cerr << "  encrypt: " << argv[0] << " encrypt <input_file> <key_file>\n";
        std::cerr << "  decrypt: " << argv[0] << " decrypt <encrypted_file> <key_file> <original_file>\n";
        return 1;
    }
    std::string mode = argv[1];
    std::string infile = argv[2];
    std::string keyfile = argv[3];

    std::ifstream fin(infile, std::ios::binary | std::ios::ate);
    if (!fin) {
        std::cerr << "Cannot open input file\n";
        return 1;
    }
    size_t file_size = fin.tellg();
    if (file_size > MAX_DATA_SIZE) {
        std::cerr << "File too large (max 1MB)\n";
        return 1;
    }
    fin.seekg(0);
    std::vector<uint8_t> data(file_size);
    fin.read((char*)data.data(), file_size);//читаем full_size байт из файла в память вектора дата возвр указ на начало внутр буфера, к char* тк реад ожид чар
    fin.close();
    std::ifstream fkey(keyfile, std::ios::binary | std::ios::ate);
    if (!fkey) {
        std::cerr << "Cannot open key file\n";
        return 1;
    }
    size_t key_size = fkey.tellg();
    if (key_size > MAX_KEY_SIZE) {
        std::cerr << "Key too long (max 256 bytes)\n";
        return 1;
    }
    fkey.seekg(0);
    std::vector<uint8_t> key(key_size);
    fkey.read((char*)key.data(), key_size);
    fkey.close();
    int shm_fd = shm_open(SHM_NAME, O_RDWR, 0666);
    if (shm_fd == -1) {
        std::cerr << "shm_open failed\n";
        return 1;
    }
    slots = (slot*)mmap(0, sizeof(slot) * NUM_SLOTS, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (slots == MAP_FAILED) {
        std::cerr << "mmap failed\n";
        return 1;
    }
    close(shm_fd);
    int slot_id = -1;
    for (int i = 0; i < NUM_SLOTS; ++i) {
        if (slots[i].status == 0) {
            slot_id = i;
            break;
        }
    }
    if (slot_id == -1) {
        std::cerr << "No free slot\n";
        munmap(slots, sizeof(slot) * NUM_SLOTS);
        return 1;
    }
    std::string sname = std::string(SEM_SERVER_BASE) + std::to_string(slot_id);
    std::string cname = std::string(SEM_CLIENT_BASE) + std::to_string(slot_id);
    sem_server = sem_open(sname.c_str(), 0);//откр сущ сем не созд новый (опис для 0)
    sem_client = sem_open(cname.c_str(), 0);
    if (sem_server == SEM_FAILED || sem_client == SEM_FAILED) {
        std::cerr << "sem_open failed\n";
        munmap(slots, sizeof(slot) * NUM_SLOTS);
        return 1;
    }
    slot& s = slots[slot_id];
    s.key_len = key_size;
    memcpy((void*)s.key, key.data(), key_size);//копирует содерж вектора key в массив s.key, дата возвр указ на данные вектора
    s.data_len = file_size;
    memcpy((void*)s.data, data.data(), file_size);
    s.status = 1;
    sem_post(sem_server);//клиент уведомл сервер о сущ данных,увел знач сем сервера
    sem_wait(sem_client);//блок потока уменьшает знач сем клиента
    if (s.status != 2) {
        std::cerr << "Server error\n";
        s.status = 0;
        sem_close(sem_server);
        sem_close(sem_client);
        munmap(slots, sizeof(slot) * NUM_SLOTS);
        return 1;
    }
    std::string outfile = (mode == "encrypt") ? infile + ".enc" : infile + ".dec";//если мод зашифровка, то к имени вх ф добав суффикс иначе
    std::ofstream fout(outfile, std::ios::binary);//выход файловый поток (созд объект класса) с именем файла в бинар режиме (без преобразования симв перевода строки)
    if (!fout) {
        std::cerr << "Cannot write output file\n";
        s.status = 0;
        sem_close(sem_server);
        sem_close(sem_client);
        munmap(slots, sizeof(slot) * NUM_SLOTS);
        return 1;
    }
    fout.write((char*)s.result, s.result_len);
    fout.close();
    s.status = 0;
    sem_close(sem_server);
    sem_close(sem_client);
    munmap(slots, sizeof(slot) * NUM_SLOTS);
    if (mode == "encrypt") {
        std::cout << "Encryption done, output: " << outfile << "\n";
        std::cout << "To decrypt, run: " << argv[0] << " decrypt " << outfile << " " << keyfile << " " << infile << "\n";
    }
    else if (mode == "decrypt") {
        if (argc < 5) {
            std::cerr << "For decrypt, provide original file as fourth argument\n";
            return 1;
        }
        std::string original = argv[4];
        if (compare_files(original, outfile)) {
            std::cout << "Decryption successful: original and decrypted files match\n";
        }
        else {
            std::cout << "Decryption failed: files differ\n";
        }
    }
    return 0;
}//дескриптор номерок,указ в ос на объект
//cd / mnt / c / rc4_project
//g++ - std = c++14 - pthread - lrt - o server server.cpp
//g++ - std = c++14 - pthread - lrt - o client client.cpp
//echo "Hello, world!" > test.txt
//echo - n "mykey123" > key.bin
//. / server
//cd / mnt / c / rc4_project
//. / client encrypt test.txt key.bin
//. / client decrypt test.txt.enc key.bin test.txt