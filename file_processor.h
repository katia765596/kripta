#ifndef FILE_PROCESSOR_H
#define FILE_PROCESSOR_H
#include "crypto_modes.h"
#include <cstddef>
#include <future>//для асинхр опер
#include <string>
class file_processor
{
public:
    file_processor(
        const crypto_modes& cryptor,
        size_t thread_count = 1
    );//принимает ссылку на объект крипто_модс (настроен на определенный режим и шифр) и кол-во потоков для параллельной обработки
    bool encrypt_file(//синх шифрование
        const std::string& input_path,
        const std::string& output_path
    ) const;
    bool decrypt_file(
        const std::string& input_path,
        const std::string& output_path
    ) const;//асинхр шифр позволяет получить рез после завершения операции 
    std::future<bool> encrypt_file_async(
        const std::string& input_path,
        const std::string& output_path
    ) const;
    std::future<bool> decrypt_file_async(
        const std::string& input_path,
        const std::string& output_path
    ) const;
    size_t get_thread_count() const;
private:
    const crypto_modes& cryptor;
    size_t thread_count;
    bool process(//фактич чтение,запись и вызов метода
        const std::string& input_path,
        const std::string& output_path,
        bool encrypt
    ) const;
};
#endif
