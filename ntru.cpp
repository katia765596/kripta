#include "ntru.h"
#include <stdexcept>
#include <vector>
namespace
{const std::vector<big_int> TRINARY = { -1, 0, 1 };}//множ коэф из которых выбираются случ полиномы в NTRU
std::pair<ntru_encrypt::public_key, ntru_encrypt::private_key>
ntru_encrypt::generate_keys(const ntru_parameters& params)
{
    if (params.N < 3 || params.p <= 1 || params.q <= params.p)
    {throw std::invalid_argument("invalid NTRU parameters");}
    const int max_attempts = 100000;//попытки генерации ключей (на случай если случ полиномы не дадут обратимых элементво)
    for (int attempt = 0; attempt < max_attempts; ++attempt)
    {
        const polynomial f =
            polynomial::random(params.N - 1, TRINARY);
        const polynomial g =
            polynomial::random(params.N - 1, TRINARY);
        try
        {
            const polynomial f_p =
                f.inverse_mod(params.p, params.N);//вычисляет обратные полиномы для f по модулю p 
            const polynomial f_q =
                f.inverse_mod(params.q, params.N);//по модулю q, все в кольце R = Z[x]/(x^N - 1)
            const polynomial h =
                f_q.mul_mod(g, params.q, params.N);//перемножает полиномы и редуцирует по модулю q и x^N-1 вычисляет открытый ключ
            return {
                public_key{params, h},
                private_key{params, f, f_p}
            };
        }
        catch (const std::runtime_error&) {}
    }
    throw std::runtime_error("NTRU key generation failed");
}
ntru_encrypt::ciphertext ntru_encrypt::encrypt(
    const public_key& pub,
    const polynomial_type& message)
{
    polynomial_type msg = message;
    msg.resize(pub.params.N, 0);//обрезает или расширяет до размера N
    msg = msg.mod(pub.params.p);//приводит коэф по модулю р
    const polynomial_type r =
        polynomial::random(pub.params.N - 1, TRINARY);//разовый ключ для шифра
    const polynomial_type rh =
        r.mul_mod(pub.h, pub.params.q, pub.params.N);//rh = r * h mod (q, x^N - 1)
    const polynomial_type e =//вычисляем шифротекст e = (p * rh + msg) mod q
        (rh * pub.params.p + msg).mod(pub.params.q);
    return { e };
}
polynomial ntru_encrypt::decrypt(
    const private_key& priv,
    const ciphertext& cipher)//возвр полином сооб
{
    const ntru_parameters& params = priv.params;
    polynomial_type a =
        priv.f.mul_mod(cipher.e, params.q, params.N);//a = f * e mod (q, x^N - 1)
    a = a.center_lift(params.q);//центрированный подъем для устранения неоднозначности при делении на p
    const polynomial_type message =
        priv.f_p.mul_mod(a, params.p, params.N);
    return message.mod(params.p);
}
polynomial ntru_encrypt::message_from_string(
    const std::string& text,
    const ntru_parameters& params)
{
    if (params.p < 256)//р велоико для кодирования байтов (каждый коэф может хранить знач до р-1
    {throw std::invalid_argument("message_from_string requires p >= 256");}
    polynomial_type::container coefficients(params.N, 0);
    for (size_t i = 0; i < text.size() && i < params.N; ++i)//копирует байты строки в первые коэф полинома (не более N cимволов)
    {coefficients[i] =static_cast<unsigned char>(text[i]);}
    return polynomial_type(coefficients);
}
std::string ntru_encrypt::message_to_string(
    const polynomial_type& poly,
    const ntru_parameters& params)
{
    if (params.p < 256)
    {throw std::invalid_argument("message_to_string requires p >= 256");}
    const auto coefficients = poly.get_coeffs();
    std::string result;
    for (size_t i = 0; i < coefficients.size() && i < params.N; ++i)
    {
        big_int value = coefficients[i] % params.p;//берет коэф по модулю p
        if (value < 0)
        {value += params.p;}
        if (value > 255)
        {throw std::runtime_error("decoded coefficient is not a byte");}
        result.push_back(
            static_cast<char>(value.convert_to<unsigned int>()));
    }
    while (!result.empty() && result.back() == '\0')
    {result.pop_back();}
    return result;
}
