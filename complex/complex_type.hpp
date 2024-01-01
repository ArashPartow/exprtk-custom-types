/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Simple Complex Type Example                                *
 * Author: Arash Partow (1999-2024)                           *
 * URL: https://www.partow.net/programming/exprtk/index.html  *
 *                                                            *
 * Copyright notice:                                          *
 * Free use of the Mathematical Expression Toolkit Library is *
 * permitted under the guidelines and in accordance with the  *
 * most current version of the MIT License.                   *
 * https://www.opensource.org/licenses/MIT                    *
 * SPDX-License-Identifier: MIT                               *
 *                                                            *
 **************************************************************
*/


#ifndef INCLUDE_COMPLEX_TYPE_HPP
#define INCLUDE_COMPLEX_TYPE_HPP


#include <cmath>
#include <complex>
#include <limits>


namespace cmplx
{
   struct complex_t
   {
      typedef std::complex<double>  value_type;

      complex_t()
      : c_(0.0)
      {}

      complex_t(const complex_t& d)
      : c_(d.c_)
      {}

      complex_t(const double& real,const double& imag)
      : c_(real,imag)
      {}

      complex_t(const value_type& v)
      : c_(v)
      {}

      complex_t& operator=(const double &d) { c_ = d; return *this; }

      template <typename T>
      complex_t& operator=(const T d) { c_ = (double)d; return *this; }

      inline complex_t& operator  =(const complex_t& r) { c_  = r.c_; return *this; }
      inline complex_t& operator +=(const complex_t& r) { c_ += r.c_; return *this; }
      inline complex_t& operator -=(const complex_t& r) { c_ -= r.c_; return *this; }
      inline complex_t& operator *=(const complex_t& r) { c_ *= r.c_; return *this; }
      inline complex_t& operator /=(const complex_t& r) { c_ /= r.c_; return *this; }

      inline complex_t& operator  =(const double r) { c_  = r; return *this; }
      inline complex_t& operator +=(const double r) { c_ += r; return *this; }
      inline complex_t& operator -=(const double r) { c_ -= r; return *this; }
      inline complex_t& operator *=(const double r) { c_ *= r; return *this; }
      inline complex_t& operator /=(const double r) { c_ /= r; return *this; }

      inline complex_t& operator++() { c_ += 1.0; return *this; }
      inline complex_t& operator--() { c_ -= 1.0; return *this; }

      inline complex_t operator++(int) { complex_t tmp(c_); c_ += 1.0; return tmp; }
      inline complex_t operator--(int) { complex_t tmp(c_); c_ -= 1.0; return tmp; }

      inline complex_t operator-() const { complex_t d; d.c_ = -c_; return d;     }

      template <typename T>
      inline operator T()    const { return static_cast<T>(c_.real()); }
      inline operator bool() const { return (c_ != 0.0);        }

      inline bool operator ==(const complex_t& r) const { return (c_ == r.c_); }
      inline bool operator !=(const complex_t& r) const { return (c_ != r.c_); }

      std::complex<double> c_;
   };

   inline complex_t operator+(const complex_t r0, const complex_t r1) { return complex_t(r0.c_ + r1.c_); }
   inline complex_t operator-(const complex_t r0, const complex_t r1) { return complex_t(r0.c_ - r1.c_); }
   inline complex_t operator*(const complex_t r0, const complex_t r1) { return complex_t(r0.c_ * r1.c_); }
   inline complex_t operator/(const complex_t r0, const complex_t r1) { return complex_t(r0.c_ / r1.c_); }


   inline bool operator< (const complex_t c0, const complex_t c1) { return std::arg(c0.c_) <  std::arg(c1.c_); }
   inline bool operator> (const complex_t c0, const complex_t c1) { return std::arg(c0.c_) >  std::arg(c1.c_); }
   inline bool operator<=(const complex_t c0, const complex_t c1) { return std::arg(c0.c_) <= std::arg(c1.c_); }
   inline bool operator>=(const complex_t c0, const complex_t c1) { return std::arg(c0.c_) >= std::arg(c1.c_); }

   #define complex_define_inequalities(Type)                                                          \
   inline complex_t operator+ (const Type& r0, const complex_t& r1) { return complex_t(r0 + r1.c_); } \
   inline complex_t operator- (const Type& r0, const complex_t& r1) { return complex_t(r0 - r1.c_); } \
   inline complex_t operator* (const Type& r0, const complex_t& r1) { return complex_t(r0 * r1.c_); } \
   inline complex_t operator/ (const Type& r0, const complex_t& r1) { return complex_t(r0 / r1.c_); } \
   inline bool operator==(const Type& r0, const complex_t& r1) { return (r0 == r1.c_);    }           \
   inline bool operator!=(const Type& r0, const complex_t& r1) { return (r0 != r1.c_);    }           \
   inline complex_t operator+ (const complex_t& r0, const Type& r1) { return complex_t(r0.c_ + r1); } \
   inline complex_t operator- (const complex_t& r0, const Type& r1) { return complex_t(r0.c_ - r1); } \
   inline complex_t operator* (const complex_t& r0, const Type& r1) { return complex_t(r0.c_ * r1); } \
   inline complex_t operator/ (const complex_t& r0, const Type& r1) { return complex_t(r0.c_ / r1); } \

   complex_define_inequalities(double)
}

namespace std
{
   using complex_t =  cmplx::complex_t;
   template <>
   class numeric_limits<complex_t>
   {
      typedef complex_t number_complex_t;

   public:

      static const bool is_specialized = true;
      static number_complex_t (min) ()        { return complex_t(std::numeric_limits<double>::min()); }
      static number_complex_t (max) ()        { return complex_t(std::numeric_limits<double>::max()); }
      static number_complex_t lowest()        { return -(max)(); }
      static number_complex_t epsilon()       { return complex_t(std::numeric_limits<double>::epsilon      ()); }
      static number_complex_t round_error()   { return complex_t(std::numeric_limits<double>::round_error  ()); }
      static number_complex_t infinity()      { return complex_t(std::numeric_limits<double>::infinity     ()); }
      static number_complex_t quiet_NaN()     { return complex_t(std::numeric_limits<double>::quiet_NaN    ()); }
      static number_complex_t signaling_NaN() { return complex_t(std::numeric_limits<double>::signaling_NaN()); }
      static number_complex_t denorm_min()    { return complex_t(std::numeric_limits<double>::denorm_min   ()); }
      static const int digits             = std::numeric_limits<double>::digits;
      static const int digits10           = std::numeric_limits<double>::digits10;
      static const int radix              = std::numeric_limits<double>::radix;
      static const int min_exponent       = std::numeric_limits<double>::min_exponent;
      static const int min_exponent10     = std::numeric_limits<double>::min_exponent10;
      static const int max_exponent       = std::numeric_limits<double>::max_exponent;
      static const int max_exponent10     = std::numeric_limits<double>::max_exponent10;
      static const bool has_infinity      = std::numeric_limits<double>::has_infinity;
      static const bool has_quiet_NaN     = std::numeric_limits<double>::has_quiet_NaN;
      static const bool has_signaling_NaN = std::numeric_limits<double>::has_signaling_NaN;
      static const bool has_denorm_loss   = std::numeric_limits<double>::has_denorm_loss;
      static const bool is_signed         = std::numeric_limits<double>::is_signed;
      static const bool is_integer        = std::numeric_limits<double>::is_integer;
      static const bool is_exact          = std::numeric_limits<double>::is_exact;
      static const bool is_iec559         = std::numeric_limits<double>::is_iec559;
      static const bool is_bounded        = std::numeric_limits<double>::is_bounded;
      static const bool is_modulo         = std::numeric_limits<double>::is_modulo;
      static const bool traps             = std::numeric_limits<double>::traps;
      static const float_denorm_style has_denorm = std::numeric_limits<double>::has_denorm;
      static const float_round_style round_style = std::numeric_limits<double>::round_style;
   };
}

namespace cmplx
{
   namespace details
   {

      namespace constant
      {
         static const complex_t e       = complex_t( 2.718281828459045235360);
         static const complex_t pi      = complex_t( 3.141592653589793238462);
         static const complex_t pi_2    = complex_t( 1.570796326794896619231);
         static const complex_t pi_4    = complex_t( 0.785398163397448309616);
         static const complex_t pi_180  = complex_t( 0.017453292519943295769);
         static const complex_t _1_pi   = complex_t( 0.318309886183790671538);
         static const complex_t _2_pi   = complex_t( 0.636619772367581343076);
         static const complex_t _180_pi = complex_t(57.295779513082320876798);
         static const complex_t log2    = complex_t( 0.693147180559945309417);
         static const complex_t sqrt2   = complex_t( 1.414213562373095048801);
      }
   }

   inline complex_t   abs(const complex_t v) { return complex_t(std::abs  (v.c_)); }
   inline complex_t  acos(const complex_t v) { return complex_t(std::acos (v.c_)); }
   inline complex_t  asin(const complex_t v) { return complex_t(std::asin (v.c_)); }
   inline complex_t  atan(const complex_t v) { return complex_t(std::atan (v.c_)); }
   inline complex_t  ceil(const complex_t v) { return complex_t(std::ceil (v.c_.real()),std::ceil (v.c_.imag())); }
   inline complex_t   cos(const complex_t v) { return complex_t(std::cos  (v.c_)); }
   inline complex_t  cosh(const complex_t v) { return complex_t(std::cosh (v.c_)); }
   inline complex_t   exp(const complex_t v) { return complex_t(std::exp  (v.c_)); }
   inline complex_t floor(const complex_t v) { return complex_t(std::floor(v.c_.real()),std::floor(v.c_.imag())); }
   inline complex_t   log(const complex_t v) { return complex_t(std::log  (v.c_)); }
   inline complex_t log10(const complex_t v) { return complex_t(std::log10(v.c_)); }
   inline complex_t  log2(const complex_t v) { return complex_t(std::log(v.c_) / details::constant::log2.c_); }
   inline complex_t   neg(const complex_t v) { return complex_t(-1.0 * v.c_); }
   inline complex_t   pos(const complex_t v) { return v;                      }
   inline complex_t   sin(const complex_t v) { return complex_t(std::sin  (v.c_)); }
   inline complex_t  sinh(const complex_t v) { return complex_t(std::sinh (v.c_)); }
   inline complex_t  sqrt(const complex_t v) { return complex_t(std::sqrt (v.c_)); }
   inline complex_t   tan(const complex_t v) { return complex_t(std::tan  (v.c_)); }
   inline complex_t  tanh(const complex_t v) { return complex_t(std::tanh (v.c_)); }
   inline complex_t   cot(const complex_t v) { return complex_t(1.0 / std::tan(v.c_)); }
   inline complex_t   sec(const complex_t v) { return complex_t(1.0 / std::cos(v.c_)); }
   inline complex_t   csc(const complex_t v) { return complex_t(1.0 / std::sin(v.c_)); }
   inline complex_t   r2d(const complex_t v) { return complex_t(v.c_ * details::constant::_180_pi.c_); }
   inline complex_t   d2r(const complex_t v) { return complex_t(v.c_ * details::constant::pi_180.c_ ); }
   inline complex_t   d2g(const complex_t v) { return complex_t(v.c_ * (20.0/9.0)); }
   inline complex_t   g2d(const complex_t v) { return complex_t(v.c_ * (9.0/20.0)); }
   inline complex_t  notl(const complex_t v) { return complex_t(v    != complex_t(0) ? complex_t(0) : complex_t(1)); }
   inline complex_t  frac(const complex_t v) { return complex_t(v.c_.real() - static_cast<long long>(v.c_.real()));  }
   inline complex_t trunc(const complex_t v) { return complex_t((double)static_cast<long long>(v.c_.real()));        }

   inline complex_t modulus(const complex_t v0, const complex_t v1) { return complex_t(fmod(v0.c_.real() , v1.c_.real()),fmod(v0.c_.imag() , v1.c_.imag())); }
   inline complex_t     pow(const complex_t v0, const complex_t v1) { return complex_t(std::pow(v0.c_,v1.c_)            ); }
   inline complex_t    logn(const complex_t v0, const complex_t v1) { return complex_t(std::log(v0.c_) / std::log(v1.c_)); }
   inline complex_t    root(const complex_t v0, const complex_t v1) { return pow(v0,complex_t(1.0) / v1);                  }
   inline complex_t   atan2(const complex_t v0, const complex_t v1) { return complex_t(std::atan2(v0.c_.real(),v0.c_.imag()),std::atan2(v1.c_.real(),v1.c_.imag())); }
   inline complex_t     max(const complex_t v0, const complex_t v1) { return complex_t(v0 > v1 ? v0.c_ : v1.c_);           }
   inline complex_t     min(const complex_t v0, const complex_t v1) { return complex_t(v0 < v1 ? v0.c_ : v1.c_);           }

   inline bool is_true  (const complex_t v) { return (v != complex_t(0)); }
   inline bool is_false (const complex_t v) { return (v == complex_t(0)); }

   inline complex_t equal(const complex_t v0x, const complex_t v1x)
   {
      const complex_t v0 = v0x;
      const complex_t v1 = v1x;
      static const complex_t epsilon = complex_t(0.0000000001);
      return (abs(v0 - v1) <= (max(complex_t(1),max(abs(v0),abs(v1))) * epsilon)) ? complex_t(1) : complex_t(0);
   }

   inline complex_t expm1(const complex_t vx)
   {
      const complex_t v = vx;
      if (abs(v) < complex_t(0.00001))
         return complex_t(v + (0.5 * v * v));
      else
         return complex_t(exp(v) - complex_t(1));
   }

   inline complex_t nequal(const complex_t v0, const complex_t v1)
   {
      static const complex_t epsilon = complex_t(0.0000000001);
      return (abs(v0 - v1) > (max(complex_t(1),max(abs(v0),abs(v1))) * epsilon)) ? complex_t(1) : complex_t(0);
   }

   inline complex_t log1p(const complex_t v)
   {
      if (v > complex_t(-1))
      {
         if (abs(v) > complex_t(0.0001))
         {
            return log(complex_t(1) + v);
         }
         else
            return (complex_t(-0.5) * v + complex_t(1)) * v;
      }
      else
         return complex_t(std::numeric_limits<double>::quiet_NaN());
   }

   inline complex_t round(const complex_t v)
   {
      return ((v < complex_t(0)) ? ceil(v - complex_t(0.5)) : floor(v + complex_t(0.5)));
   }

   inline complex_t roundn(const complex_t v0, const complex_t v1)
   {
      const complex_t p10 = pow(complex_t(10),trunc(v1));
      if (v0 < complex_t(0))
         return complex_t(ceil ((v0 * p10) - complex_t(0.5)) / p10);
      else
         return complex_t(floor((v0 * p10) + complex_t(0.5)) / p10);
   }

   inline complex_t hypot(const complex_t v0, const complex_t v1)
   {
      return sqrt((v0 * v0) + (v1 * v1));
   }

   inline complex_t shr(const complex_t v0, const complex_t v1)
   {
      return v0 * (complex_t(1) / pow(complex_t(2),trunc(v1)));
   }

   inline complex_t shl(const complex_t v0, const complex_t v1)
   {
      return v0 * pow(complex_t(2),trunc(v1));
   }

   inline complex_t sgn(const complex_t v)
   {
           if (v > complex_t(0)) return complex_t(+1);
      else if (v < complex_t(0)) return complex_t(-1);
      else               return complex_t( 0);
   }

   inline complex_t nand(const complex_t v0, const complex_t& v1)
   {
      return (is_false(v0) || is_false(v1)) ? complex_t(1) : complex_t(0);
   }

   inline complex_t nor(const complex_t v0, const complex_t& v1)
   {
      return (is_false(v0) && is_false(v1)) ? complex_t(1) : complex_t(0);
   }

   inline complex_t xnor(const complex_t v0, const complex_t& v1)
   {
      const bool v0_true = is_true(v0);
      const bool v1_true = is_true(v1);
      if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
         return complex_t(1);
      else
         return complex_t(0);
   }

   inline complex_t erf(complex_t /*v*/)
   {
      // Note: Implementation for erf of a complex number is required.
      // http://ab-initio.mit.edu/Faddeeva.hh
      return complex_t(0);
   }

   inline complex_t erfc(complex_t v)
   {
      return complex_t(1) - erf(v);
   }
}

#endif
