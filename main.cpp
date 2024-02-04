//#include "bigint.hpp"

#include <assert.h>
#include <exception>
#include <iostream>
#include <stdlib.h>

#include "gmp-6.3.0/gmpxx.h"

//using ix = BigInt;


class GmpWrapper {
public:
    // Constructors
    GmpWrapper() {
        mpz_init(value_);
    }

    GmpWrapper(const char* str) {
        mpz_init_set_str(value_, str, 10);
    }

    GmpWrapper(const GmpWrapper& other) {
        mpz_init_set(value_, other.value_);
    }

    GmpWrapper(int intValue) {
        mpz_init_set_si(value_, intValue);
    }

    // Destructor
    ~GmpWrapper() {
        mpz_clear(value_);
    }

    // Assignment operator
    GmpWrapper& operator=(const GmpWrapper& other) {
        if (this != &other) {
            mpz_set(value_, other.value_);
        }
        return *this;
    }

    // Arithmetic operators
    GmpWrapper operator+(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_add(result.value_, value_, other.value_);
        return result;
    }

    GmpWrapper operator-(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_sub(result.value_, value_, other.value_);
        return result;
    }

    GmpWrapper operator*(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_mul(result.value_, value_, other.value_);
        return result;
    }

    GmpWrapper operator%(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_mod(result.value_, value_, other.value_);
        return result;
    }

    GmpWrapper operator/(const GmpWrapper& other) const {
        if (mpz_sgn(other.value_) == 0) {
            // Handle division by zero
            throw std::invalid_argument("Division by zero");
        }

        GmpWrapper result;
        mpz_tdiv_q(result.value_, value_, other.value_);
        return result;
    }

        // Bitwise AND operator
    GmpWrapper operator&(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_and(result.value_, value_, other.value_);
        return result;
    }

    // Bitwise OR operator
    GmpWrapper operator|(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_ior(result.value_, value_, other.value_);
        return result;
    }

    // Bitwise XOR operator
    GmpWrapper operator^(const GmpWrapper& other) const {
        GmpWrapper result;
        mpz_xor(result.value_, value_, other.value_);
        return result;
    }

    // Bitwise NOT operator
    GmpWrapper operator~() const {
        GmpWrapper result;
        mpz_com(result.value_, value_);
        return result;
    }

    GmpWrapper pow(unsigned int exp) const {
        GmpWrapper result;
        mpz_pow_ui(result.value_, value_, exp);
        return result;
    }
  
    unsigned char get_bit(size_t bitIndex) const {
        return mpz_tstbit(value_, bitIndex);
    }

    size_t count_trailing_zeros() const {
        return static_cast<size_t>(mpz_scan1(value_, 0));
    }

    size_t bitlength() const {
        return static_cast<size_t>(mpz_sizeinbase(value_, 2));
    }

    void write() const {
        char* str = mpz_get_str(nullptr, 10, value_);
        std::cout << str;
        free(str);
    }

        // Equality operator
    bool operator==(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) == 0;
    }

    // Inequality operator
    bool operator!=(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) != 0;
    }

    // Less than operator
    bool operator<(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) < 0;
    }

    // Less than or equal to operator
    bool operator<=(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) <= 0;
    }

    // Greater than operator
    bool operator>(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) > 0;
    }

    // Greater than or equal to operator
    bool operator>=(const GmpWrapper& other) const {
        return mpz_cmp(value_, other.value_) >= 0;
    }

        // Addition with int
    GmpWrapper operator+(int intValue) const {
        GmpWrapper n(intValue);
        return this->operator+(n);
    }

    // Subtraction with int
    GmpWrapper operator-(int intValue) const {
        GmpWrapper n(intValue);
        return this->operator-(n);
    }

    GmpWrapper operator-() const {
        GmpWrapper result;
        mpz_neg(result.value_, value_);
        return result;
    }

    // Multiplication with int
    GmpWrapper operator*(int intValue) const {
        GmpWrapper result;
        mpz_mul_si(result.value_, value_, intValue);
        return result;
    }

    // Division with int
    GmpWrapper operator/(int intValue) const {
        if (intValue == 0) {
            // Handle division by zero
            throw std::invalid_argument("Division by zero");
        }

        GmpWrapper n(intValue);
        return this->operator/(n);
    }

    // Output operator
    friend std::ostream& operator<<(std::ostream& os, const GmpWrapper& gmp) {
        return os << mpz_get_str(nullptr, 10, gmp.value_);
    }

private:
    mpz_t value_;
};

using ix = GmpWrapper;


ix
ix_abs(const ix& _n)
{
  return _n < 0 ? _n * -1 : _n;
}

[[gnu::pure]] std::pair<ix, ix>
divmod(const ix& _na, const ix& _nb) noexcept
{
  if (_na > 0 && _nb > 0)
  {
    return {_na / _nb, _na - (_na / _nb) * _nb};
  }
  if (_na > 0 && _nb < 0)
  {
    const ix abs_dist = ix_abs(_nb - _na) - 1;
    const ix div      = (abs_dist / _nb) * -1;
    const ix rem      = ix_abs(_na - (div * _nb * -1));

    return {div * -1, rem * -1};
  }
  if (_na < 0 && _nb > 0)
  {
    const ix abs_dist = ix_abs(_nb - _na) - 1;
    const ix div      = (abs_dist / _nb);
    const ix rem      = ix_abs((_na * -1) - (div * _nb));

    return {div * -1, rem};
  }
  if (_na < 0 && _nb < 0)
  {
    return {_na / _nb, ((-_na) - ((-_na) / (-_nb)) * (-_nb)) * -1};
  }

  return {0, 0};
}

[[gnu::pure]] std::pair<ix, ix>
extended_gcd(const ix& aa, const ix& bb) noexcept
{
  ix lastremainder = ix_abs(aa);
  ix remainder     = ix_abs(bb);

  ix x, lasty = 0;
  ix y, lastx = 1;

  ix quotient;

  while (remainder != 0)
  {
    ix oldremainder = remainder;

    const auto dm = divmod(lastremainder, remainder);

    quotient  = dm.first;
    remainder = dm.second;

    lastremainder = oldremainder;

    const ix tmp_x = x;
    x              = lastx - quotient * x;
    lastx          = tmp_x;

    const ix tmp_y = y;
    y              = lasty - quotient * y;
    lasty          = tmp_y;
  }

  return {lastremainder, lastx * (aa < 0 ? -1 : 1)};
}

[[gnu::pure]] ix
mod(const ix& _na, const ix& _nb) noexcept
{
  ix result = _na % _nb;

  if (result < ix{0}) [[unlikely]]
  {
    return result + _nb;
  }

  return result;
}

[[gnu::pure]] ix
modinv(const ix& a, const ix& m) noexcept
{
  const auto gx = extended_gcd(a, m);
  return mod(gx.second, m);
}

struct crv_p
{
  ix x{}, y{};

  void
  print() const
  {
    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << "]\n";
  }
};

[[gnu::pure]] inline crv_p
point_add(const crv_p& _p1, const crv_p& _p2, const ix& _p /* global modulo */) noexcept
{
  const ix s = mod((_p2.y - _p1.y) * modinv(_p2.x - _p1.x, _p), _p);
  ix x3      = mod((s.pow(2) - _p1.x - _p2.x), _p);
  ix y3      = mod(s * (_p1.x - x3) - _p1.y, _p);

  return {x3, y3};
}

[[gnu::pure]] inline crv_p
point_double(const crv_p& _p1, const ix& _p) noexcept
{
  const ix x2 = _p1.x;
  ix s        = mod(((ix{3} * _p1.x.pow(2)) * modinv(ix{3} * _p1.y, _p)), _p);
  ix x3       = mod((s.pow(2) - _p1.x - x2), _p);
  ix y3       = mod((s * (_p1.x - x3) - _p1.y), _p);

  return {x3, y3};
}

[[gnu::pure]] inline std::size_t
bits_to_represent(ix _num) noexcept
{
  return _num.bitlength() - _num.count_trailing_zeros();
}

[[gnu::pure]] inline bool
is_bit_set(const ix& _num, std::size_t _bid) noexcept
{
  return bool(_num.get_bit(_bid));
}

/* TODO: stupid naive and slow */
/* possibly look at:

https://doi.org/10.1016/j.jksuci.2019.07.013
https://link.springer.com/chapter/10.1007/978-3-540-73074-3_15

*/

crv_p
double_and_add(const crv_p& _p1, const ix& _num, const ix& _p)
{
  const std::size_t bits = bits_to_represent(_num);
  const crv_p _p2        = _p1;
  crv_p _p3              = _p1;

  for (std::size_t i = 2; i != bits + 1; ++i)
  {
    if (is_bit_set(_num, bits - i))
    {
      _p3 = point_double(_p3, _p);
      _p3 = point_add(_p2, _p3, _p);
    }
    else
    {
      _p3 = point_double(_p3, _p);
    }
  }
  return _p3;
}

/*
  BIG DISCLAIMER: UNFINISHED AND NOT TESTED!!!

  TODO:
    - add better scaler multiply
    - use actually usable BigInt libary that has reasonable performance

*/

int
main()
{
  ix mod_global = "115792089237316195423570985008687907853269984665640564039457584007908834671663";

  crv_p G = {"55066263022277343669578718895168534326250603453777594175500187360389116729240",
             "32670510020758816978083085130507043184471273380659243275938904335757337482424"};

  ix privKeyA = "40505654708211189456746820883201845994248137211058198699828051064905928553035";
  ix privKeyB = "83862260130769358743610306176715755043868098730045613807339143668249321773381";

  crv_p pubKeyA = double_and_add(G, privKeyA, mod_global);
  pubKeyA.print();

  crv_p pubKeyB = double_and_add(G, privKeyB, mod_global);
  pubKeyB.print();

  crv_p shared_secretA = double_and_add(pubKeyB, privKeyA, mod_global);
  crv_p shared_secretB = double_and_add(pubKeyA, privKeyB, mod_global);

  std::cout << "shared secrets: \n";

  shared_secretA.print();
  shared_secretB.print();

  return 0;
}