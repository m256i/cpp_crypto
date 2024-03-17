// #include "bigint.hpp"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <vector>
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
  assert(a != 0);

  const auto gx = extended_gcd(a, m);
  return mod(gx.second, m);
}

struct crv_p
{
  ix x{}, y{};

  ix debug_value{};

  void
  print() const
  {

    debug_value.write();

    std::cout << " : ";

    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << "]\n";
  }

  bool
  operator==(const crv_p& other) const
  {
    return (x == other.x) && (y == other.y);
  }

  bool
  operator!=(const crv_p& other) const
  {
    return !(this->operator==(other));
  }
};

struct jcbn_crv_p
{
  ix x{}, y{}, z{};

  ix debug_value{};

  void
  print() const
  {
    debug_value.write();

    std::cout << " : ";

    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << ", ";
    z.write();
    std::cout << "]\n";
  }

  bool
  operator==(const jcbn_crv_p& other) const
  {
    return (x == other.x) && (y == other.y) && (z == other.z);
  }

  bool
  operator!=(const jcbn_crv_p& other) const
  {
    return !(this->operator==(other));
  }
};

static const crv_p a_identity_element      = {0, 0};
static const jcbn_crv_p j_identity_element = {0, 0, 1};

jcbn_crv_p
to_jacobian(const crv_p& _ws_point)
{
  return {_ws_point.x, _ws_point.y, ix{1}, _ws_point.debug_value};
}

crv_p
from_jacobian(const jcbn_crv_p& _jcbn, const ix& _p)
{
  if (_jcbn == j_identity_element)
  {
    return {0, 0};
  }

  ix inv = modinv(_jcbn.z, _p);
  return {mod(_jcbn.x * inv.pow(2), _p), mod(_jcbn.y * inv.pow(3), _p), _jcbn.debug_value};
}

/* NOTE when using curves where a != 0 this needs to be changed */

jcbn_crv_p
point_double(const jcbn_crv_p& _p1, const ix& _p)
{
  if (_p1 == jcbn_crv_p{0, 0, 1})
  {
    std::cout << "recieved O in doubling!\n";
    return {0, 0, 1};
  }

  const ix a = ((_p1.x * 4) * _p1.y.pow(2));
  const ix b = (_p1.x.pow(2) * 3);

  std::cout << "------------- point double -------------\n";

  std::cout << "  a: ";
  a.write();
  std::cout << "\n";
  std::cout << "  b: ";
  b.write();
  std::cout << "\n\n";

  jcbn_crv_p out;

  out.x = ((a * -2) + b.pow(2)) % _p;

  std::cout << "  x: ";
  out.x.write();
  std::cout << "\n";

  out.y = ((_p1.y.pow(4) * -8) + b * (a - out.x)) % _p;

  std::cout << "  y: ";
  out.y.write();
  std::cout << "\n";

  out.z = (_p1.y * _p1.z * 2) % _p;
  if (out.z == 0)
  {
    out.z = 1;
  }

  std::cout << "  z: ";
  out.z.write();
  std::cout << "\n";

  out.debug_value = (_p1.debug_value * 2) % _p;

  return out;
}

/* TODO: fix */
jcbn_crv_p
point_add(const jcbn_crv_p& _p1, const jcbn_crv_p& _p2, const ix& _p)
{
  if (_p1 == _p2)
  {
    std::cout << "recieved same point in addition!\n";
    return point_double(_p1, _p);
  }

  if (_p1 == jcbn_crv_p{0, 0, 1})
  {
    std::cout << "recieved O p1 in addition!\n";
    return _p2;
  }

  if (_p2 == jcbn_crv_p{0, 0, 1})
  {
    std::cout << "recieved O p2 in addition!\n";
    return _p1;
  }

  assert(_p1 != j_identity_element && _p2 != j_identity_element);

  const ix U1 = _p1.x * _p2.z.pow(2);
  const ix U2 = _p2.x * _p1.z.pow(2);
  const ix S1 = _p1.y * _p2.z.pow(3);
  const ix S2 = _p2.y * _p1.z.pow(3);

  jcbn_crv_p out;

  const ix H = (U2 - U1);
  const ix R = (S2 - S1);

  out.x = (R.pow(2) - H.pow(3) - 2 * U1 * H.pow(2)) % _p;
  out.y = (R * (U1 * H.pow(2) - out.x) - S1 * H.pow(3)) % _p;
  out.z = (H * _p1.z * _p2.z) % _p;

  if (out.z == 0)
  {
    out.z = 1;
  }

  out.debug_value = (_p1.debug_value + _p2.debug_value) % _p;

  return out;
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

[[gnu::pure]] inline crv_p
point_add(const crv_p& _p1, const crv_p& _p2, const ix& _p) noexcept
{
  if (_p1 == _p2)
  {
    return point_double(_p1, _p);
  }

  const ix s = mod((_p2.y - _p1.y) * modinv(_p2.x - _p1.x, _p), _p);
  ix x3      = mod((s.pow(2) - _p1.x - _p2.x), _p);
  ix y3      = mod(s * (_p1.x - x3) - _p1.y, _p);
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

SIKE all of these papers had errors and typos, i wrote this shit myself :/
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

  for (std::size_t i = 1; i != bits; ++i)
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

static constexpr std::size_t window_size = 4;

std::vector<jcbn_crv_p>
precompute(const jcbn_crv_p& Q, const ix& _p)
{
  const std::size_t count = (std::size_t)std::pow(2, window_size);
  std::vector<jcbn_crv_p> out;
  out.reserve(count);

  out.push_back({0, 0, 1});
  out.push_back(Q);

  jcbn_crv_p next{Q}; // 2P

  for (std::size_t i = 0; i != count; ++i)
  {
    std::cout << "adding: \n  ";
    Q.print();
    std::cout << "  ";
    next.print();

    next = point_add(Q, next, _p); // ++P
    out.push_back(next);

    std::cout << "  result: ";
    next.print();
  }

  std::printf("size of precomp table: %lu\n", out.size());
  for (std::size_t i = 0; i != out.size(); ++i)
  {
    std::cout << i << " G: ";
    from_jacobian(out[i], _p).print();
    out[i].print();
  }

  // precomp now {O, 1P, 2P, 3P, ..., (2^w-1)P} -> size 16 for w = 4

  return out;
}

jcbn_crv_p
windowed_scalar_mul(const std::vector<jcbn_crv_p>& _precomp, const ix& _num, const ix& _p)
{
  jcbn_crv_p Q{0, 0, 1};

  int current_q = 0;

  std::size_t m = bits_to_represent(_num) / window_size;

  while ((m * window_size) < bits_to_represent(_num))
  {
    ++m;
  }

  // std::printf("bits_to_represent(_num): %lu\n", bits_to_represent(_num));
  // std::printf("w: %lu\n", window_size);
  // std::printf("m: %lu\n", m);
  // std::printf("precomp size: %lu\n", _precomp.size());
  // _num.writeb();
  // std::puts("");

  for (std::size_t i = 0; i != m; ++i)
  {

    std::cout << current_q << "Q ->";
    from_jacobian(Q, _p).print();

    const std::size_t nbits = _num.get_bits(((signed)m - (signed)i - 1) * window_size, window_size);

    if (nbits > 0)
    {

      Q = point_add(Q, _precomp[nbits], _p);
      current_q += nbits;

      std::cout << "Q added " << nbits << " -> " << current_q << "Q ";
      from_jacobian(Q, _p).print();

      if (i < m - 1)
      {
        for (auto j = 0ul; j != window_size; ++j)
        {
          Q = point_double(Q, _p);
          current_q *= 2;

          std::cout << "Q doubled "
                    << " -> " << current_q << "Q ";
          from_jacobian(Q, _p).print();
        }
      }
    }

    // std::printf("(m-i-1) %lu\n", ((signed)m - (signed)i - 1) * window_size);
    // std::cout << std::bitset<4>(nbits) << "\n";
  }

  std::puts("");

  return Q;
}

int
main()
{
  // ix mod_global = "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";

  // crv_p G = {"0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
  //            "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"};

  // ix privKeyA = "40505654708211189456746820883201845994248137211058198699828051064905928553035";
  // ix privKeyB = "83862260130769358743610306176715755043868098730045613807339143668249321773381";

  ix mod_global{17};
  jcbn_crv_p G = {15, 13, 1, 1};

  ix privKeyA = 5;
  ix privKeyB = 9;

  const auto G_precomp = precompute(G, mod_global);

  jcbn_crv_p pubKeyA = windowed_scalar_mul(G_precomp, privKeyA, mod_global);
  jcbn_crv_p pubKeyB = windowed_scalar_mul(G_precomp, privKeyB, mod_global);

  std::cout << "--------- pubkeys ---------\n";

  from_jacobian(pubKeyA, mod_global).print();
  pubKeyA.print();
  from_jacobian(pubKeyB, mod_global).print();
  pubKeyB.print();

  std::cout << "---------------------------\n";

  jcbn_crv_p pubKeyAJ = pubKeyA;
  jcbn_crv_p pubKeyBJ = pubKeyB;

  std::cout << "jacobian: \n";
  {
    perf_ _("jacobian");
    jcbn_crv_p shared_secretAJ = double_and_add(pubKeyBJ, privKeyA, mod_global);
    jcbn_crv_p shared_secretBJ = double_and_add(pubKeyAJ, privKeyB, mod_global);
    from_jacobian(shared_secretAJ, mod_global).print();
    from_jacobian(shared_secretBJ, mod_global).print();
  }

  const auto precomp = precompute(pubKeyBJ, mod_global);
  // const auto precomp2 = precompute(pubKeyBJ, mod_global);

  std::cout << "jacobian windowed: \n";
  {
    perf_ _("jacobian windowed");

    jcbn_crv_p shared_secretAJ = windowed_scalar_mul(precomp, privKeyA, mod_global);
    from_jacobian(shared_secretAJ, mod_global).print();

    // jcbn_crv_p shared_secretBJ = windowed_scalar_mul(precomp2, pubKeyAJ, privKeyB, mod_global);
    // from_jacobian(shared_secretBJ, mod_global).print();
  }

  // std::cout << "affine: \n";
  // {
  //   perf_ _("affine");
  //   crv_p shared_secretA = double_and_add(pubKeyB, privKeyA, mod_global);
  //   crv_p shared_secretB = double_and_add(pubKeyA, privKeyB, mod_global);
  //   shared_secretA.print();
  //   shared_secretB.print();
  // }

  return 0;
}

// here be dragons!