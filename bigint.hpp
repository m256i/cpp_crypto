/*
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>
© https://github.com/983/bigint
*/

#ifndef NUM_HPP_INCLUDED
#define NUM_HPP_INCLUDED

#include "bigint.h"

#include <stdlib.h>

struct BigInt
{
  bigint data[1];

  BigInt() { bigint_init(data); }

  BigInt(int b)
  {
    bigint_init(data);
    bigint_from_int(data, b);
  }

  BigInt(const char *s, int base = 10)
  {
    bigint_init(data);
    bigint_from_str_base(data, s, base);
  }

  BigInt(const BigInt &b)
  {
    bigint_init(data);
    bigint_cpy(data, b.data);
  }

  BigInt &
  operator=(const BigInt &b)
  {
    bigint_cpy(data, b.data);
    return *this;
  }

  ~BigInt() { bigint_free(data); }

  static void
  div_mod(BigInt &quotient, BigInt &remainder, const BigInt &biginterator, const BigInt &denominator)
  {
    bigint_div_mod(quotient.data, remainder.data, biginterator.data, denominator.data);
  }

  void
  write(char *dst, int *n_dst, int dst_base = 10, int zero_terminate = 1) const
  {
    bigint_write_base(dst, n_dst, data, dst_base, zero_terminate);
  }

  template <class STREAM>
  STREAM &
  write(STREAM &s, int dst_base = 10) const
  {
    int n     = bigint_write_size(data, dst_base);
    char *buf = (char *)malloc(n);
    write(buf, &n, dst_base);
    s << buf;
    free(buf);
    return s;
  }

  BigInt &
  operator<<=(int shift)
  {
    bigint_shift_left(data, data, shift);
    return *this;
  }

  BigInt &
  operator>>=(int shift)
  {
    bigint_shift_right(data, data, shift);
    return *this;
  }

  BigInt &
  operator+=(const BigInt &b)
  {
    bigint_add(data, data, b.data);
    return *this;
  }

  BigInt &
  operator-=(const BigInt &b)
  {
    bigint_sub(data, data, b.data);
    return *this;
  }

  BigInt &
  operator*=(const BigInt &b)
  {
    bigint_mul(data, data, b.data);
    return *this;
  }

  BigInt &
  operator/=(const BigInt &b)
  {
    bigint_div(data, data, b.data);
    return *this;
  }

  BigInt &
  operator%=(const BigInt &b)
  {
    bigint_mod(data, data, b.data);
    return *this;
  }

  BigInt &
  operator++()
  {
    bigint_add_word(data, data, 1);
    return *this;
  }

  BigInt &
  operator--()
  {
    bigint_sub_word(data, data, 1);
    return *this;
  }

  BigInt &
  set_bit(int bit_index)
  {
    bigint_set_bit(data, bit_index);
    return *this;
  }

  BigInt &
  clr_bit(int bit_index)
  {
    bigint_clr_bit(data, bit_index);
    return *this;
  }

  bigint_word
  get_bit(int bit_index) const
  {
    return bigint_get_bit(data, bit_index);
  }

  int
  bitlength() const
  {
    return bigint_bitlength(data);
  }

  int
  count_trailing_zeros() const
  {
    return bigint_count_trailing_zeros(data);
  }

  bool
  is_probable_prime(int n_tests, bigint_rand_func rand_func) const
  {
    return bigint_is_probable_prime(data, n_tests, rand_func);
  }

  BigInt
  sqrt() const
  {
    BigInt b;
    bigint_sqrt(b.data, data);
    return b;
  }

  BigInt
  pow(bigint_word exponent) const
  {
    BigInt b;
    bigint_pow_word(b.data, data, exponent);
    return b;
  }

  static BigInt
  gcd(const BigInt &a, const BigInt &b)
  {
    BigInt c;
    bigint_gcd(c.data, a.data, b.data);
    return c;
  }

  static BigInt
  rand_bits(int n_bits, bigint_rand_func rand_func)
  {
    BigInt b;
    bigint_rand_bits(b.data, n_bits, rand_func);
    return b;
  }

  static BigInt
  rand_inclusive(const BigInt &n, bigint_rand_func rand_func)
  {
    BigInt b;
    bigint_rand_inclusive(b.data, n.data, rand_func);
    return b;
  }

  static BigInt
  rand_exclusive(const BigInt &n, bigint_rand_func rand_func)
  {
    BigInt b;
    bigint_rand_exclusive(b.data, n.data, rand_func);
    return b;
  }
};

inline BigInt
operator-(const BigInt &a)
{
  BigInt b(a);
  b.data->neg = !b.data->neg;
  return b;
}

inline BigInt
operator+(const BigInt &a, const BigInt &b)
{
  return BigInt(a) += b;
}

inline BigInt
operator-(const BigInt &a, const BigInt &b)
{
  return BigInt(a) -= b;
}

inline BigInt
operator*(const BigInt &a, const BigInt &b)
{
  return BigInt(a) *= b;
}

inline BigInt
operator/(const BigInt &a, const BigInt &b)
{
  return BigInt(a) /= b;
}

inline BigInt
operator%(const BigInt &a, const BigInt &b)
{
  return BigInt(a) %= b;
}

inline BigInt
operator<<(const BigInt &a, int shift)
{
  return BigInt(a) <<= shift;
}

inline BigInt
operator>>(const BigInt &a, int shift)
{
  return BigInt(a) >>= shift;
}

inline bool
operator==(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) == 0;
}
inline bool
operator!=(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) != 0;
}
inline bool
operator<=(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) <= 0;
}
inline bool
operator>=(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) >= 0;
}
inline bool
operator<(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) < 0;
}
inline bool
operator>(const BigInt &a, const BigInt &b)
{
  return bigint_cmp(a.data, b.data) > 0;
}

#endif
