#include "file_processor.h"
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <vector>
file_processor::file_processor(const crypto_modes& cryptor_value,size_t thread_count_value)
    : cryptor(cryptor_value),thread_count(thread_count_value)
{
    if (thread_count == 0) {throw std::invalid_argument("invalid thread count");}
}
bool file_processor::process(const std::string& input_path, const std::string& output_path, bool encrypt
) const
{
    try
    {
        std::ifstream input(input_path, std::ios::binary);
        if (!input) {return false;}
        byte_array source(
            (std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>()
        );
        input.close();
        byte_array result;
        if (encrypt)
        {
            result = cryptor.encrypt_parallel(source,thread_count);
        }
        else
        {
            result = cryptor.decrypt_parallel(source,thread_count);
        }
        std::ofstream output(output_path,std::ios::binary | std::ios::trunc);//в бин режиме, транк- если файл сущ то его содержимое будет усечено(перезаписано)
        if (!output)
        {
            return false;
        }
        output.write(
            reinterpret_cast<const char*>(result.data()),//возвр указ на начало массива байтов, котрый приводится к конст чар
            static_cast<std::streamsize>(result.size())
        );
        if (!output)
        {
            return false;
        }
        output.close();
        return true;
    }
    catch (...)
    {
        return false;
    }
}
bool file_processor::encrypt_file(//синхронное шифр файла
    const std::string& input_path,
    const std::string& output_path
) const
{
    return process(
        input_path,
        output_path,
        true
    );
}
bool file_processor::decrypt_file(
    const std::string& input_path,
    const std::string& output_path
) const
{
    return process(
        input_path,
        output_path,
        false
    );
}
std::future<bool> file_processor::encrypt_file_async(
    const std::string& input_path,
    const std::string& output_path
) const
{
    return std::async(
        std::launch::async,
        [this, input_path, output_path]()
        {
            return encrypt_file(
                input_path,
                output_path
            );
        }
    );
}
std::future<bool> file_processor::decrypt_file_async(
    const std::string& input_path,
    const std::string& output_path
) const
{
    return std::async(
        std::launch::async,//гарантирует асинхронное выполнение в отдельном потоке
        [this, input_path, output_path]()
        {
            return decrypt_file(
                input_path,
                output_path
            );
        }
    );
}
size_t file_processor::get_thread_count() const
{
    return thread_count;
}
