#ifndef BIG_INT_H
#define BIG_INT_H
#include <boost/multiprecision/cpp_int.hpp>
using big_int = boost::multiprecision::cpp_int;
big_int random_bigint(const big_int& min, const big_int& max);
big_int mod_pow(const big_int& base, const big_int& exp, const big_int& mod);
big_int mod_inverse(const big_int& a, const big_int& mod);
bool is_prime(const big_int& n, int certainty = 25);//число итераций вер-го теста чем выше, тем меньше ошибка
big_int generate_prime(int bits);
big_int center_lift(const big_int& x, const big_int& mod);
#endif
