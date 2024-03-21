// #include "bigint.hpp"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <vector>
#include <bitset>
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
    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << "]";

    std::cout << " (";
    debug_value.write();
    std::cout << ")\n";
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

  void
  print() const
  {
    std::cout << "[";
    x.write();
    std::cout << ", ";
    y.write();
    std::cout << ", ";
    z.write();
    std::cout << "]";

    // std::cout << " (";
    // debug_value.write();
    // std::cout << ")\n";
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
static const jcbn_crv_p j_identity_element = {1, 1, 0};

jcbn_crv_p
to_jacobian(const crv_p& _ws_point)
{
  return {_ws_point.x, _ws_point.y, ix{1}};
}

crv_p
from_jacobian(const jcbn_crv_p& _jcbn, const ix& _p)
{
  if (_jcbn == j_identity_element ||_jcbn.z == 0)
  {
    return {0, 0};
  }

  ix inv = modinv(_jcbn.z, _p);
  return {mod(_jcbn.x * inv.pow(2), _p), mod(_jcbn.y * inv.pow(3), _p)};
}

/* NOTE when using curves where a != 0 this needs to be changed */

jcbn_crv_p
point_double(const jcbn_crv_p& _p1, const ix& _p)
{
  if (_p1.y == 0) [[unlikely]] {
    // std::cout << "y = 0 in point doubling!\n";
    return j_identity_element;
  }

  const ix a = (_p1.x * 4) * _p1.y.pow(2) % _p;
  const ix b = _p1.x.pow(2) * 3 % _p /* + a * _p1.z.pow(4) */;

  jcbn_crv_p out;

  out.x = (b.pow(2) - 2 * a) % _p;
  out.y = ((_p1.y.pow(4) * -8) + b * (a - out.x)) % _p;
  out.z = (_p1.y * _p1.z * 2) % _p;

  //out.debug_value = (_p1.debug_value * 2) /* % _p*/;

  return out;
}

/* TODO: fix */
jcbn_crv_p
point_add(const jcbn_crv_p& _p1, const jcbn_crv_p& _p2, const ix& _p)
{
  if (_p1 == j_identity_element) 
  {
    return {_p2.x, _p2.y, _p2.z};
  }
  else if (_p2 == j_identity_element)
  {
    return {_p1.x, _p1.y, _p1.z};
  }

  const ix U1 = _p1.x * _p2.z.pow(2);
  const ix U2 = _p2.x * _p1.z.pow(2);
  const ix S1 = _p1.y * _p2.z.pow(3);
  const ix S2 = _p2.y * _p1.z.pow(3);

  if (U1 == U2) [[unlikely]]
  {
    if (S1 != S2) /*very*/[[unlikely]] {
      return j_identity_element;
    }
    else {
      return point_double(_p1, _p);
    }
  }

  jcbn_crv_p out;

  const ix H = U2 - U1;
  const ix R = S2 - S1;

  out.x = (R.pow(2) - H.pow(3) - 2 * U1 * H.pow(2)) % _p;
  out.y = (R * (U1 * H.pow(2) - out.x) - S1 * H.pow(3)) % _p;
  out.z = (H * _p1.z * _p2.z) % _p;
  
  if (out.z == 0) [[unlikely]]
  {
    out = j_identity_element;
  }

  return out;
}

[[gnu::pure]] inline std::size_t
bits_to_represent(ix _num) noexcept
{
  return _num.bitlength(); 
}

[[gnu::pure]] inline bool
is_bit_set(const ix& _num, std::size_t _bid) noexcept
{
  return bool(_num.get_bit(_bid));
}

/* possibly look at:

https://doi.org/10.1016/j.jksuci.2019.07.013
https://link.springer.com/chapter/10.1007/978-3-540-73074-3_15

SIKE all of these papers had errors and typos, i wrote this shit myself :/
*/

static constexpr std::size_t window_size = 4;

std::vector<jcbn_crv_p>
precompute(const jcbn_crv_p& Q, const ix& _p)
{
  const std::size_t count = (std::size_t)std::pow(2, window_size);
  std::vector<jcbn_crv_p> out;
  out.reserve(count);
  
  out.push_back(j_identity_element);
  out.push_back(Q);

  jcbn_crv_p next{Q}; // 1P

  for (std::size_t i = 2; i != count; ++i)
  {
    next = point_add(Q, next, _p); // ++P
    out.push_back(next);
  }
  // precomp now {O, 1P, 2P, 3P, ..., (2^w-1)P} -> size 16 for w = 4

  return out;
}

jcbn_crv_p
windowed_scalar_mul(const std::vector<jcbn_crv_p>& _precomp, const ix& _num, const ix& _p)
{
  jcbn_crv_p Q{j_identity_element};
  std::size_t m = bits_to_represent(_num) / window_size;

  while ((m * window_size) < bits_to_represent(_num))
  {
    ++m;
  }

  for (std::size_t i = 0; i != m; ++i)
  {
    for (auto j = 0ul; j != window_size; ++j)
    {
      Q = point_double(Q, _p);
    }

    const std::size_t start_idx = ((signed)m - (signed)i - 1) * window_size;
    const std::size_t nbits = _num.get_bits(start_idx, window_size);

    if (nbits > 0) [[likely]]
    {
      Q = point_add(Q, _precomp[nbits], _p);
    }
  }

  return Q;
}

int
main()
{
  ix mod_global = "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";
  crv_p G = {"0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
             "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 1};
  // ix privKeyA = "0x598D635BD02C77CC3020CFFD744D4D75D190C41E726D16C2FE2F5A1F06AC324B";
  // ix privKeyB = "0xb9685b6ee0405eb5389c9b9d29404357eec208f05471b21e58dad170371f9945";

  ix privKeyA = "0x598D635BD02C77CC3020CFFD744D4D75D190C41E726D16C2FE2F5A1F06AC324B";
  ix privKeyB = "0xb9685b6ee0405eb5389c94897d04357eec208f05471b21e58dad170371f9945";
  const auto G_precomp = precompute(to_jacobian(G), mod_global);

  //ix mod_global{17};
  //jcbn_crv_p G = {15, 13, 1, 1};
  //ix privKeyA = 32121984718;
  //ix privKeyB = 67222;
  //const auto G_precomp = precompute(G, mod_global);

  jcbn_crv_p pubKeyA = windowed_scalar_mul(G_precomp, privKeyA, mod_global);
  jcbn_crv_p pubKeyB = windowed_scalar_mul(G_precomp, privKeyB, mod_global);


  std::cout << "--------- pubkeys ---------\n";

  //pubKeyA.print();
  from_jacobian(pubKeyA, mod_global).print();
  from_jacobian(pubKeyB, mod_global).print();
  // pubKeyB.print();
  std::cout << "---------------------------\n";

  const jcbn_crv_p pubKeyAJ = pubKeyA;
  const jcbn_crv_p pubKeyBJ = pubKeyB;

  const auto precomp = precompute(pubKeyBJ, mod_global);
  const auto precomp2 = precompute(pubKeyAJ, mod_global);
  // const auto precomp2 = precompute(pubKeyBJ, mod_global);

 jcbn_crv_p shared_secretAJ, shared_secretBJ;

  std::cout << "jacobian windowed: \n";
  {
    perf_ _("jacobian windowed");

    shared_secretAJ = windowed_scalar_mul(precomp, privKeyA, mod_global);
    shared_secretBJ = windowed_scalar_mul(precomp2, privKeyB, mod_global);
  }
  
  from_jacobian(shared_secretAJ, mod_global).print();
  from_jacobian(shared_secretBJ, mod_global).print();

  // std::cout << "jacobian: \n";
  // {
  //   perf_ _("jacobian");
  //   jcbn_crv_p shared_secretAJ = double_and_add(pubKeyBJ, privKeyA, mod_global);
  //   jcbn_crv_p shared_secretBJ = double_and_add(pubKeyAJ, privKeyB, mod_global);
  //   from_jacobian(shared_secretAJ, mod_global).print();
  //   from_jacobian(shared_secretBJ, mod_global).print();
  // }

  return 0;
}

// here be dragons!
