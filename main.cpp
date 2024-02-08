// #include "bigint.hpp"

#include <assert.h>
#include <exception>
#include <iostream>
#include <stdlib.h>
#include <chrono>

#include "crypto.h"

using namespace blue_crypto;

using ix = GmpWrapper;

struct perf_
{
  perf_(std::string_view _name) : name(_name) { start = std::chrono::high_resolution_clock::now(); }

  ~perf_()
  {
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << name << ": " << duration.count() << " microseconds" << std::endl;
  }

  std::string_view name;
  decltype(std::chrono::high_resolution_clock::now()) start;
};

ix
ix_abs(const ix& _n)
{
  return _n < 0 ? _n * -1 : _n;
}

[[gnu::pure]] ix
mod(const ix& _na, const ix& _nb) noexcept
{
  ix result = _na % _nb;
  return result;
}

[[gnu::pure]] std::pair<ix, ix>
divmod(const ix& _na, const ix& _nb) noexcept
{
  return {_na / _nb, mod(_na, _nb)};
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

struct jcbn_crv_p
{
  ix x{}, y{}, z{};

  void
  print() const
  {
    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << ", ";
    z.write();
    std::cout << "]\n";
  }
};

jcbn_crv_p
to_jacobian(const crv_p& _ws_point)
{
  return {_ws_point.x, _ws_point.y, ix{1}};
}

crv_p
from_jacobian(const jcbn_crv_p& _jcbn, const ix& _p)
{
  ix inv = modinv(_jcbn.z, _p);
  return {mod(_jcbn.x * inv.pow(2), _p), mod(_jcbn.y * inv.pow(3), _p)};
}

/* NOTE when using curves where a != 0 this needs to be changed */

jcbn_crv_p
point_double(const jcbn_crv_p& _p1, const ix& _p)
{
  const ix a = mod((_p1.x * 4) * _p1.y.pow(2), _p);
  const ix b = mod(_p1.x.pow(2) * 3, _p);

  jcbn_crv_p out;

  out.x = mod((a * -2) + b.pow(2), _p);
  out.y = mod((_p1.y.pow(4) * -8) + b * (a - out.x), _p);
  out.z = mod(_p1.y * _p1.z * 2, _p);

  return out;
}

/* TODO: fix */

jcbn_crv_p
point_add(const jcbn_crv_p& _p1, const jcbn_crv_p& _p2, const ix& _p)
{
  const ix a = mod(_p1.x * _p1.z.pow(2), _p);
  const ix b = mod(_p2.x * _p1.z.pow(2), _p);
  const ix c = mod(_p1.y * _p2.z.pow(3), _p);
  const ix d = mod(_p2.y * _p1.z.pow(3), _p);

  const ix e = mod(b - a, _p);
  const ix f = mod(d - c, _p);

  jcbn_crv_p out;

  out.x = (-e.pow(3) - 2 * a * e.pow(2) + f) % _p;
  out.y = -c * e.pow(3) + f * (a * e.pow(2) - out.x);
  out.z = (_p1.z * _p2.z * e) % _p;

  return out;
}

[[gnu::pure]] inline crv_p
point_add(const crv_p& _p1, const crv_p& _p2, const ix& _p) noexcept
{
  const ix s = mod((_p2.y - _p1.y) * modinv(_p2.x - _p1.x, _p), _p);
  ix x3      = mod((s.pow(2) - _p1.x - _p2.x), _p);
  ix y3      = mod(s * (_p1.x - x3) - _p1.y, _p);
  return {x3, y3};
}

[[gnu::pure]] inline crv_p
point_double(const crv_p& _p1, const ix& _p) noexcept
{
  const ix& x2 = _p1.x;

  ix s  = mod(((ix{3} * _p1.x.pow(2)) * modinv(ix{2} * _p1.y, _p)), _p);
  ix x3 = mod((s.pow(2) - _p1.x - x2), _p);
  ix y3 = mod((s * (_p1.x - x3) - _p1.y), _p);

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

jcbn_crv_p
double_and_add(const jcbn_crv_p& _p1, const ix& _num, const ix& _p)
{
  const std::size_t bits = bits_to_represent(_num);
  const jcbn_crv_p _p2   = _p1;
  jcbn_crv_p _p3         = _p1;

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
  ix mod_global = "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";

  crv_p G = {"0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
             "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"};

  ix privKeyA = "40505654708211189456746820883201845994248137211058198699828051064905928553035";
  ix privKeyB = "83862260130769358743610306176715755043868098730045613807339143668249321773381";

  auto p2G = point_double(G, mod_global);

  std::cout << "singular double of G: ";
  {
    perf_ _("affine");

    // point_double(point_double(G, mod_global), mod_global).print();
    point_add(G, p2G, mod_global).print();
  }

  std::cout << "jacobian: ";
  {
    perf_ _("jacobian");
    // from_jacobian(point_double(point_double(to_jacobian(G), mod_global), mod_global), mod_global).print();
    from_jacobian(point_add(to_jacobian(G), to_jacobian(p2G), mod_global), mod_global).print();
  }

  // crv_p pubKeyA = double_and_add(G, privKeyA, mod_global);
  // pubKeyA.print();
  // crv_p pubKeyB = double_and_add(G, privKeyB, mod_global);
  // pubKeyB.print();
  // jcbn_crv_p pubKeyAJ = to_jacobian(pubKeyA);
  // jcbn_crv_p pubKeyBJ = to_jacobian(pubKeyB);
  // jcbn_crv_p shared_secretAJ = double_and_add(pubKeyBJ, privKeyA, mod_global);
  // jcbn_crv_p shared_secretBJ = double_and_add(pubKeyAJ, privKeyB, mod_global);
  // crv_p shared_secretA = double_and_add(pubKeyB, privKeyA, mod_global);
  // crv_p shared_secretB = double_and_add(pubKeyA, privKeyB, mod_global);
  // std::cout << "shared secrets: \n";
  // shared_secretA.print();
  // shared_secretB.print();
  // from_jacobian(shared_secretAJ, mod_global).print();
  // from_jacobian(shared_secretBJ, mod_global).print();

  return 0;
}