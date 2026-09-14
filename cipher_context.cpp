#include "cipher_context.h"
#include "modes.h"
#include "paddings.h"
#include <fstream>
#include <thread>
#include <future>
#include <stdexcept>
#include <vector>
#include <cstring>
cipher_context::cipher_context(i_symmetric_cipher* algo,
    cipher_mode mode,
    padding_mode pad,
    const std::vector<uint8_t>& iv,
    const std::initializer_list<uint8_t>& extra)
    : cipher(algo), iv(iv), extra_params(extra) {
    if (!cipher) throw std::invalid_argument("null cipher");
    block_size = cipher->block_size();
    switch (mode) {
    case cipher_mode::ecb: this->mode.reset(new ecb_mode()); break;//заменяет тек умн указ новым объектом автом удаляя старый
    case cipher_mode::cbc: this->mode.reset(new cbc_mode()); break;
    case cipher_mode::pcbc: this->mode.reset(new pcbc_mode()); break;
    case cipher_mode::cfb: this->mode.reset(new cfb_mode()); break;
    case cipher_mode::ofb: this->mode.reset(new ofb_mode()); break;
    case cipher_mode::ctr: this->mode.reset(new ctr_mode()); break;
    case cipher_mode::random_delta: this->mode.reset(new random_delta_mode()); break;
    default: throw std::invalid_argument("unknown mode");
    }
    switch (pad) {
    case padding_mode::zeros: this->padding.reset(new zeros_padding()); break;
    case padding_mode::ansi_x923: this->padding.reset(new ansi_x923_padding()); break;
    case padding_mode::pkcs7: this->padding.reset(new pkcs7_padding()); break;
    case padding_mode::iso10126: this->padding.reset(new iso10126_padding()); break;
    default: throw std::invalid_argument("unknown padding");
    }
    this->mode->init(iv, extra_params);//иниц режим передавая ему iv и доп параметры
}
cipher_context::~cipher_context() {}//шифер не удалится (контекст не владеет алг)
void cipher_context::encrypt(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, int num_threads) {
    process_blocks(input, output, true, num_threads);
}
void cipher_context::decrypt(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, int num_threads) {
    process_blocks(input, output, false, num_threads);
}
std::future<void> cipher_context::encrypt_async(const std::string& input_path, const std::string& output_path, int num_threads) {
    return std::async(std::launch::async, [this, input_path, output_path, num_threads]() {
        process_file_sync(input_path, output_path, true, num_threads);
        });//асунк запускает переданную лямбу в отдельном потоке флаг лаун гарантирует асинхр выполнение лямба захватывает вис, пути к файлам и число потоков, внутри синхр методо возвр фьюча который позволяет ожид завершегия операции
}
std::future<void> cipher_context::decrypt_async(const std::string& input_path, const std::string& output_path, int num_threads) {
    return std::async(std::launch::async, [this, input_path, output_path, num_threads]() {
        process_file_sync(input_path, output_path, false, num_threads);
        });
}
void cipher_context::process_blocks(const std::vector<uint8_t>& input, std::vector<uint8_t>& output, bool encrypt, int num_threads) {
    if (input.empty()) { output.clear(); return; }
    if (num_threads <= 0) num_threads = 1;
    mode->reset();
    mode->init(iv, extra_params);
    std::vector<uint8_t> data = input;
    bool use_padding = mode->needs_padding();
    if (encrypt) {//паддинг для блочных режимов при шифр и дешифр
        if (use_padding) {
            data = padding->add_padding(data, block_size);
        }
    }
    else {
        if (use_padding) {
            if (data.size() % block_size != 0) throw std::invalid_argument("ciphertext size not multiple of block");
        }
    }
    size_t data_len = data.size();
    size_t full_blocks = data_len / block_size;
    size_t remainder = data_len % block_size;//неполный блок для блочных с паддингом =0
    bool parallel = encrypt ? mode->can_parallel_encrypt() : mode->can_parallel_decrypt();
    if (parallel && num_threads > 1 && full_blocks > 0) {//если да, то запускаем многопоточную обработку
        std::vector<std::thread> threads;
        std::vector<std::vector<uint8_t>> results(full_blocks);//один вектор - один блок
        size_t blocks_per_thread = (full_blocks + num_threads - 1) / num_threads;//скок блоков на каждый поток
        for (int t = 0; t < num_threads; ++t) {
            threads.emplace_back([this, &data, &results, full_blocks, blocks_per_thread, t, num_threads, encrypt]() {//создает  новый пток прямо в векторе (доб в конец) создает лямбду которая будет выполнена в новом потоке, знает свой номер t и диапозон блоков которые должна обработать
                size_t start = t * blocks_per_thread;//номер потока на кол-во потоков на один блок
                size_t end = std::min(start + blocks_per_thread, full_blocks);//вычислили индексы блоков для тек потока
                for (size_t i = start; i < end; ++i) {//цикл по блокам в потоке
                    std::vector<uint8_t> block(block_size);
                    std::memcpy(block.data(), data.data() + i * block_size, block_size);//извл 1 блок из вх данных и помещаем его в отдельный вектор
                    results[i] = mode->process_block(block, encrypt, cipher);//вызов вирт метода у объекта мод
                }
                });//закрывает лямбду и передаем в констр потока 
        }
        for (auto& th : threads) th.join();//идем по векторам потоков джоин блок тек поток до завершении этого потока
        output.clear();
        output.reserve(full_blocks * block_size + remainder);
        for (size_t i = 0; i < full_blocks; ++i) {//эл в конец вектора
            output.insert(output.end(), results[i].begin(), results[i].end());
        }
    }
    else {
        output.clear();//если параллелить незя то посл в осн потоке обраб каждый блок
        output.reserve(full_blocks * block_size + remainder);
        for (size_t i = 0; i < full_blocks; ++i) {
            std::vector<uint8_t> block(block_size);
            std::memcpy(block.data(), data.data() + i * block_size, block_size);
            std::vector<uint8_t> processed = mode->process_block(block, encrypt, cipher);
            output.insert(output.end(), processed.begin(), processed.end());
        }
    }
    if (remainder > 0 && !use_padding) {
        std::vector<uint8_t> last_block(block_size, 0);
        std::memcpy(last_block.data(), data.data() + full_blocks * block_size, remainder);//копируем остаток байтов в начало блока (остальное 0)
        std::vector<uint8_t> processed = mode->process_block(last_block, encrypt, cipher);
        output.insert(output.end(), processed.begin(), processed.begin() + remainder);
    }
    if (!encrypt && use_padding) {
        output = padding->remove_padding(output, block_size);
    }
}
void cipher_context::process_file_sync(const std::string& input_path, const std::string& output_path, bool encrypt, int num_threads) {//синхр обработка файлов метод вызывается из асинхронных оберток
    std::ifstream fin(input_path, std::ios::binary);
    if (!fin) throw std::runtime_error("cannot open input file");
    fin.seekg(0, std::ios::end);//указ от конца фвйла
    size_t size = fin.tellg();//возвр тек поз указ (размер файла в байтах)
    fin.seekg(0, std::ios::beg);
    std::vector<uint8_t> data(size);
    if (size > 0) fin.read(reinterpret_cast<char*>(data.data()), size);//читаем сайз байт из файла в буфер вектора
    fin.close();
    std::vector<uint8_t> result;
    process_blocks(data, result, encrypt, num_threads);
    std::ofstream fout(output_path, std::ios::binary);
    if (!fout) throw std::runtime_error("cannot open output file");
    fout.write(reinterpret_cast<const char*>(result.data()), result.size());
    fout.close();
}