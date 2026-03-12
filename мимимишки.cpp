#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cstdint>
#include <thread>
#include <future>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include <random>
#include <locale>
namespace fs = std::filesystem; 
using bytearr = std::vector<uint8_t>;
enum class shifrmode { ECB, CBC, PCBC };
enum class paddingmod { Zeros, ANSI_X923 };
class Isymmshifr {
public:
    virtual ~Isymmshifr() = default;
    virtual size_t bl_size() const = 0;
    virtual void encryptb(const uint8_t* input, uint8_t* output, const bytearr& key) = 0;
    virtual void decryptb(const uint8_t* input, uint8_t* output, const bytearr& key) = 0;
};
class Imode {
public:
    virtual ~Imode() = default;
    virtual void encrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) = 0;
    virtual void decrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) = 0;
};
class Ipad {
public:
    virtual ~Ipad() = default;
    virtual void addp(bytearr& data, size_t bsize) = 0;
    virtual void rempad(bytearr& data, size_t bsize) = 0;
};
class Zerospad : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t padlenn = bsize - (data.size() % bsize);
        if (padlenn == 0) padlenn = bsize;
        data.insert(data.end(), padlenn, 0);
    }
    void rempad(bytearr& data, size_t) override {
        while (!data.empty() && data.back() == 0) data.pop_back();
    }
};
class ANSIX923Padding : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t padlenn = bsize - (data.size() % bsize);
        if (padlenn == 0) padlenn = bsize;
        data.insert(data.end(), padlenn - 1, 0);
        data.push_back(static_cast<uint8_t>(padlenn));
    }
    void rempad(bytearr& data, size_t) override {
        if (data.empty()) return;
        uint8_t padlenn = data.back();
        if (padlenn > data.size()) throw std::runtime_error("Invalid padding");
        data.erase(data.end() - padlenn, data.end());
    }
};
class ECBMode : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr&, size_t thr) override {
        size_t bsize = alg->bl_size();
        output = data;
        size_t allb = data.size() / bsize;
        std::vector<std::thread> thrs;
        for (size_t t = 0; t < thr; ++t) {
            thrs.emplace_back([=, &output]() {
                for (size_t i = t; i < allb; i += thr)
                    alg->encryptb(&data[i * bsize], &output[i * bsize], key);
                });
        }
        for (auto& th : thrs) th.join();
    }
    void decrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr&, size_t thr) override {
        size_t bsize = alg->bl_size();
        output = data;
        size_t allb = data.size() / bsize;
        std::vector<std::thread> thrs;
        for (size_t t = 0; t < thr; ++t) {
            thrs.emplace_back([=, &output]() {
                for (size_t i = t; i < allb; i += thr)
                    alg->decryptb(&data[i * bsize], &output[i * bsize], key);
                });
        }
        for (auto& th : thrs) th.join();
    }
};
class CBCMode : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t bsize = alg->bl_size();
        if (iv.size() != bsize) throw std::runtime_error("Invalid for CBC");
        output.resize(data.size());
        bytearr prbl = iv;
        for (size_t i = 0; i < data.size(); i += bsize) {
            bytearr block(data.begin() + i, data.begin() + i + bsize);
            for (size_t j = 0; j < bsize; ++j) block[j] ^= prbl[j];
            alg->encryptb(block.data(), &output[i], key);
            prbl.assign(output.begin() + i, output.begin() + i + bsize);
        }
    }
    void decrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t bsize = alg->bl_size();
        if (iv.size() != bsize) throw std::runtime_error("Invalid for CBC");
        output.resize(data.size());
        bytearr prbl = iv;
        for (size_t i = 0; i < data.size(); i += bsize) {
            bytearr tempb(data.begin() + i, data.begin() + i + bsize);
            alg->decryptb(&data[i], &output[i], key);
            for (size_t j = 0; j < bsize; ++j) output[i + j] ^= prbl[j];
            prbl = tempb;
        }
    }
};
class PCBCMode : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t bsize = alg->bl_size();
        if (iv.size() != bsize) throw std::runtime_error("Invalid for PCBC");
        output.resize(data.size());
        bytearr prbl = iv;
        for (size_t i = 0; i < data.size(); i += bsize) {
            bytearr block(data.begin() + i, data.begin() + i + bsize);
            for (size_t j = 0; j < bsize; ++j) block[j] ^= prbl[j];
            alg->encryptb(block.data(), &output[i], key);
            for (size_t j = 0; j < bsize; ++j) prbl[j] = block[j] ^ output[i + j];
        }
    }
    void decrypt(const bytearr& data, bytearr& output, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t bsize = alg->bl_size();
        if (iv.size() != bsize) throw std::runtime_error("Invalid IV for PCBC");
        output.resize(data.size());
        bytearr prbl = iv;
        for (size_t i = 0; i < data.size(); i += bsize) {
            bytearr cipherbl(data.begin() + i, data.begin() + i + bsize);
            alg->decryptb(&data[i], &output[i], key);
            for (size_t j = 0; j < bsize; ++j) output[i + j] ^= prbl[j];
            for (size_t j = 0; j < bsize; ++j) prbl[j] = output[i + j] ^ cipherbl[j];
        }
    }
};
class Symm {
public:
    Symm(std::shared_ptr<Isymmshifr> alg,
        shifrmode encmode,
        paddingmod padmode,
        const bytearr& iv = {})
        : algorithm_(alg), iv_(iv) {
        padding_ = createpadding(padmode);
        mode_ = createmode(encmode, iv_);
    }
    void encrypt(const bytearr& input, bytearr& output, const bytearr& key, size_t thr = 1) {
        bytearr data = input;
        padding_->addp(data, algorithm_->bl_size());
        mode_->encrypt(data, output, algorithm_.get(), key, iv_, thr);
    }
    void decrypt(const bytearr& input, bytearr& output, const bytearr& key, size_t thr = 1) {
        mode_->decrypt(input, output, algorithm_.get(), key, iv_, thr);
        padding_->rempad(output, algorithm_->bl_size());
    }
    std::future<void> encryptfasy(const std::string& infile, const std::string& outfile,
        const bytearr& key, size_t thr = 1) {
        return std::async(std::launch::async, [=]() {
            std::ifstream in(infile, std::ios::binary);
            std::ofstream out(outfile, std::ios::binary);
            if (!in || !out) throw std::runtime_error("File open failed");
            bytearr buffer((std::istreambuf_iterator<char>(in)), {});
            bytearr encrypted;
            encrypt(buffer, encrypted, key, thr);
            out.write(reinterpret_cast<char*>(encrypted.data()), encrypted.size());
            });
    }
    std::future<void> decryptfasy(const std::string& infile, const std::string& outfile,
        const bytearr& key, size_t thr = 1) {
        return std::async(std::launch::async, [=]() {
            std::ifstream in(infile, std::ios::binary);
            std::ofstream out(outfile, std::ios::binary);
            if (!in || !out) throw std::runtime_error("File open failed");
            bytearr buffer((std::istreambuf_iterator<char>(in)), {});
            bytearr decrypted;
            decrypt(buffer, decrypted, key, thr);
            out.write(reinterpret_cast<char*>(decrypted.data()), decrypted.size());
            });
    }
private:
    std::shared_ptr<Isymmshifr> algorithm_;
    std::shared_ptr<Imode> mode_;
    std::shared_ptr<Ipad> padding_;
    bytearr iv_;
    std::shared_ptr<Imode> createmode(shifrmode mode, const bytearr& iv) {
        switch (mode) {
        case shifrmode::ECB: return std::make_shared<ECBMode>();
        case shifrmode::CBC: return std::make_shared<CBCMode>();
        case shifrmode::PCBC: return std::make_shared<PCBCMode>();
        default: throw std::runtime_error("Unsupported encryption mode");
        }
    }
    std::shared_ptr<Ipad> createpadding(paddingmod padMode) {
        switch (padMode) {
        case paddingmod::Zeros: return std::make_shared<Zerospad>();
        case paddingmod::ANSI_X923: return std::make_shared<ANSIX923Padding>();
        default: throw std::runtime_error("Unsupported padding mode");
        }
    }
};
class XORalg : public Isymmshifr {
public:
    size_t bl_size() const override { return 8; }
    void encryptb(const uint8_t* input, uint8_t* output, const bytearr& key) override {
        for (size_t i = 0; i < bl_size(); ++i)
            output[i] = input[i] ^ key[i % key.size()];
    }
    void decryptb(const uint8_t* input, uint8_t* output, const bytearr& key) override {
        encryptb(input, output, key);
    }
};
bytearr randkey(size_t length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, 255);
    bytearr key(length);
    for (auto& b : key) b = static_cast<uint8_t>(dis(gen));
    return key;
}
void unittest() {
    std::cout << "Запуск unit-тестов...\n";
    auto algo = std::make_shared<XORalg>();
    bytearr key = randkey(8);
    bytearr iv(8, 0);
    Zerospad zp;
    bytearr data = { 1,2,3 };
    zp.addp(data, 8);
    zp.rempad(data, 8);
    if (data != bytearr({ 1,2,3 })) throw std::runtime_error("ZerosPadding test failed");
    ANSIX923Padding ap;
    data = { 4,5,6 };
    ap.addp(data, 8);
    ap.rempad(data, 8);
    if (data != bytearr({ 4,5,6 })) throw std::runtime_error("ANSIX923Padding test failed");
    bytearr block(8, 0xAA);
    bytearr out(8);
    algo->encryptb(block.data(), out.data(), key);
    bytearr decrypted(8);
    algo->decryptb(out.data(), decrypted.data(), key);
    if (block != decrypted) throw std::runtime_error("Block encrypt/decrypt test failed");
    bytearr datamode = { 1,2,3,4,5,6,7,8 };
    bytearr encrypted;
    bytearr decrypteddata;
    ECBMode ecb; CBCMode cbc; PCBCMode pcbc;
    std::vector<Imode*> modes = { &ecb, &cbc, &pcbc };
    for (auto m : modes) {
        m->encrypt(datamode, encrypted, algo.get(), key, iv, 2);
        m->decrypt(encrypted, decrypteddata, algo.get(), key, iv, 2);
        if (decrypteddata != datamode) throw std::runtime_error("Mode test failed");
    }
    fs::create_directories("test_input");
    std::ofstream("test_input/test1.txt") << "test file content";
    Symm ctx(algo, shifrmode::ECB, paddingmod::Zeros);
    ctx.encryptfasy("test_input/test1.txt", "test_enc.enc", key, 2).get();
    ctx.decryptfasy("test_enc.enc", "test_dec.txt", key, 2).get();
    std::ifstream orig("test_input/test1.txt"), dec("test_dec.txt");
    std::vector<char> o((std::istreambuf_iterator<char>(orig)), {});
    std::vector<char> d((std::istreambuf_iterator<char>(dec)), {});
    if (o != d) throw std::runtime_error("File test failed");
    fs::remove_all("test_input");
    fs::remove("test_enc.enc");
    fs::remove("test_dec.txt");
    std::cout << "Все unit-тесты пройдены!\n";
}
int main() {
    setlocale(LC_ALL, "Russian");
    try {
        auto algo = std::make_shared<XORalg>();
        bytearr key = randkey(8);
        std::cout << "Сгенерирован ключ: ";
        for (auto b : key) std::cout << static_cast<int>(b) << " ";
        std::cout << "\n";
        bytearr iv(8, 0);
        int modeinput;
        std::cout << "Выберите режим шифрования (0-ECB, 1-CBC, 2-PCBC): ";
        std::cin >> modeinput;
        shifrmode mode = static_cast<shifrmode>(modeinput);
        int padinput;
        std::cout << "Выберите паддинг (0-Zeros, 1-ANSI X.923): ";
        std::cin >> padinput;
        paddingmod pad = static_cast<paddingmod>(padinput);
        Symm ctx(algo, mode, pad, iv);
        std::string fpath, be, outf;
        std::cout << "Введите путь к папке с файлами для шифрования: ";
        std::cin >> fpath;
        std::cout << "Папка для зашифрованных файлов: ";
        std::cin >> be;
        std::cout << "Папка для расшифрованных файлов: ";
        std::cin >> outf;
        fs::create_directories(be);
        fs::create_directories(outf);
        std::vector<std::future<void>> futures;
        for (auto& entry : fs::directory_iterator(fpath)) {
            if (!entry.is_regular_file()) continue;
            std::string inf = entry.path().string();
            std::string encryptedf = (fs::path(be) / entry.path().filename()).string();
            std::string decryptedf = (fs::path(outf) / entry.path().filename()).string();
            futures.push_back(std::async(std::launch::async, [&, inf, encryptedf, decryptedf]() {
                ctx.encryptfasy(inf, encryptedf, key).get();
                ctx.decryptfasy(encryptedf, decryptedf, key).get();
                std::cout << "Файл обработан: " << entry.path().filename() << std::endl;
                }));
        }
        for (auto& fut : futures) fut.get();
        unittest();
    }
    catch (const std::exception& ex) {
        std::cerr << "Ошибка: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
