/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Simple Real Type Example                                   *
 * Authors: Arash Partow                                      *
 * URL: http://www.partow.net/programming/exprtk/index.html   *
 *                                                            *
 * Copyright notice:                                          *
 * Free use of the Mathematical Expression Toolkit Library is *
 * permitted under the guidelines and in accordance with the  *
 * most current version of the MIT License.                   *
 * http://www.opensource.org/licenses/MIT                     *
 *                                                            *
 **************************************************************
*/


#ifndef INCLUDE_REAL_TYPE_HPP
#define INCLUDE_REAL_TYPE_HPP


#include <cmath>
#include <algorithm>
#include <limits>


namespace real
{
   struct type
   {
      type()
      : d_(0.0)
      {}

      type(const type& d)
      : d_(d.d_)
      {}

      type(const double& d)
      : d_(d)
      {}

      template <typename T>
      explicit type(const T& d)
      : d_((double)d)
      {}

      type& operator=(const double &d) { d_ = d; return *this; }

      template <typename T>
      type& operator=(const T d) { d_ = (double)d; return *this; }

      inline type& operator  =(const type& r) { d_  = r.d_; return *this; }
      inline type& operator +=(const type& r) { d_ += r.d_; return *this; }
      inline type& operator -=(const type& r) { d_ -= r.d_; return *this; }
      inline type& operator *=(const type& r) { d_ *= r.d_; return *this; }
      inline type& operator /=(const type& r) { d_ /= r.d_; return *this; }

      inline type& operator  =(const double r) { d_  = r; return *this; }
      inline type& operator +=(const double r) { d_ += r; return *this; }
      inline type& operator -=(const double r) { d_ -= r; return *this; }
      inline type& operator *=(const double r) { d_ *= r; return *this; }
      inline type& operator /=(const double r) { d_ /= r; return *this; }

      inline type& operator++() { d_ += 1.0; return *this; }
      inline type& operator--() { d_ -= 1.0; return *this; }

      inline type operator++(int) { type tmp(d_); d_ += 1.0; return tmp; }
      inline type operator--(int) { type tmp(d_); d_ -= 1.0; return tmp; }

      inline type operator-() const { type d; d.d_ = -d_; return d;     }

      template <typename T>
      inline operator T()    const { return static_cast<T>(d_); }
      inline operator bool() const { return (d_ != 0.0);        }

      inline bool operator ==(const type& r) const { return (d_ == r.d_); }
      inline bool operator !=(const type& r) const { return (d_ != r.d_); }

      double d_;
   };

   inline type operator+(const type r0, const type r1) { return type(r0.d_ + r1.d_); }
   inline type operator-(const type r0, const type r1) { return type(r0.d_ - r1.d_); }
   inline type operator*(const type r0, const type r1) { return type(r0.d_ * r1.d_); }
   inline type operator/(const type r0, const type r1) { return type(r0.d_ / r1.d_); }

   inline bool operator< (const type r0, const type r1) { return (r0.d_ <  r1.d_); }
   inline bool operator> (const type r0, const type r1) { return (r0.d_ >  r1.d_); }
   inline bool operator<=(const type r0, const type r1) { return (r0.d_ <= r1.d_); }
   inline bool operator>=(const type r0, const type r1) { return (r0.d_ >= r1.d_); }

   #define real_define_inequalities(Type)                                                      \
   inline type operator+ (const Type/*&*/ r0, const type/*&*/ r1) { return type(r0 + r1.d_); } \
   inline type operator- (const Type/*&*/ r0, const type/*&*/ r1) { return type(r0 - r1.d_); } \
   inline type operator* (const Type/*&*/ r0, const type/*&*/ r1) { return type(r0 * r1.d_); } \
   inline type operator/ (const Type/*&*/ r0, const type/*&*/ r1) { return type(r0 / r1.d_); } \
   inline bool operator< (const Type/*&*/ r0, const type/*&*/ r1) { return (r0 <  r1.d_);    } \
   inline bool operator> (const Type/*&*/ r0, const type/*&*/ r1) { return (r0 >  r1.d_);    } \
   inline bool operator<=(const Type/*&*/ r0, const type/*&*/ r1) { return (r0 <= r1.d_);    } \
   inline bool operator>=(const Type/*&*/ r0, const type/*&*/ r1) { return (r0 >= r1.d_);    } \
   inline bool operator==(const Type/*&*/ r0, const type/*&*/ r1) { return (r0 == r1.d_);    } \
   inline bool operator!=(const Type/*&*/ r0, const type/*&*/ r1) { return (r0 != r1.d_);    } \
   inline type operator+ (const type/*&*/ r0, const Type/*&*/ r1) { return type(r0.d_ + r1); } \
   inline type operator- (const type/*&*/ r0, const Type/*&*/ r1) { return type(r0.d_ - r1); } \
   inline type operator* (const type/*&*/ r0, const Type/*&*/ r1) { return type(r0.d_ * r1); } \
   inline type operator/ (const type/*&*/ r0, const Type/*&*/ r1) { return type(r0.d_ / r1); } \
   inline bool operator< (const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ <  r1);    } \
   inline bool operator> (const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ >  r1);    } \
   inline bool operator<=(const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ <= r1);    } \
   inline bool operator>=(const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ >= r1);    } \
   inline bool operator==(const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ == r1);    } \
   inline bool operator!=(const type/*&*/ r0, const Type/*&*/ r1) { return (r0.d_ != r1);    } \

   real_define_inequalities(double      )
   real_define_inequalities(float       )
   real_define_inequalities(int         )
   real_define_inequalities(long long   )
   real_define_inequalities(unsigned int)
   real_define_inequalities(unsigned long long)
   real_define_inequalities(unsigned long int )
}

namespace std
{
   template <>
   class numeric_limits<real::type>
   {
      typedef real::type number_type;

   public:

      static const bool is_specialized = true;
      static number_type (min) ()        { return std::numeric_limits<double>::min(); }
      static number_type (max) ()        { return std::numeric_limits<double>::max(); }
      static number_type lowest()        { return -(max)(); }
      static number_type epsilon()       { return std::numeric_limits<double>::epsilon();       }
      static number_type round_error()   { return std::numeric_limits<double>::round_error();   }
      static number_type infinity()      { return std::numeric_limits<double>::infinity();      }
      static number_type quiet_NaN()     { return std::numeric_limits<double>::quiet_NaN();     }
      static number_type signaling_NaN() { return std::numeric_limits<double>::signaling_NaN(); }
      static number_type denorm_min()    { return std::numeric_limits<double>::denorm_min();    }
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

namespace real
{
   namespace details
   {
      namespace constant
      {
         static const type e       =  2.718281828459045235360;
         static const type pi      =  3.141592653589793238462;
         static const type pi_2    =  1.570796326794896619231;
         static const type pi_4    =  0.785398163397448309616;
         static const type pi_180  =  0.017453292519943295769;
         static const type _1_pi   =  0.318309886183790671538;
         static const type _2_pi   =  0.636619772367581343076;
         static const type _180_pi = 57.295779513082320876798;
         static const type log2    =  0.693147180559945309417;
         static const type sqrt2   =  1.414213562373095048801;
      }
   }

   inline type   abs(const type v) { return std::abs  (v.d_); }
   inline type  acos(const type v) { return std::acos (v.d_); }
   inline type  asin(const type v) { return std::asin (v.d_); }
   inline type  atan(const type v) { return std::atan (v.d_); }
   inline type  ceil(const type v) { return std::ceil (v.d_); }
   inline type   cos(const type v) { return std::cos  (v.d_); }
   inline type  cosh(const type v) { return std::cosh (v.d_); }
   inline type   exp(const type v) { return std::exp  (v.d_); }
   inline type floor(const type v) { return std::floor(v.d_); }
   inline type   log(const type v) { return std::log  (v.d_); }
   inline type log10(const type v) { return std::log10(v.d_); }
   inline type  log2(const type v) { return std::log(v.d_) / type(details::constant::log2); }
   inline type   neg(const type v) { return type(-1.0) * v;   }
   inline type   pos(const type v) { return v;                }
   inline type   sin(const type v) { return std::sin  (v.d_); }
   inline type  sinh(const type v) { return std::sinh (v.d_); }
   inline type  sqrt(const type v) { return std::sqrt (v.d_); }
   inline type   tan(const type v) { return std::tan  (v.d_); }
   inline type  tanh(const type v) { return std::tanh (v.d_); }
   inline type   cot(const type v) { return type(1) / std::tan(v.d_); }
   inline type   sec(const type v) { return type(1) / std::cos(v.d_); }
   inline type   csc(const type v) { return type(1) / std::sin(v.d_); }
   inline type   r2d(const type v) { return (v.d_ * type(details::constant::_180_pi)); }
   inline type   d2r(const type v) { return (v.d_ * type(details::constant::pi_180));  }
   inline type   d2g(const type v) { return (v.d_ * type(20.0/9.0)); }
   inline type   g2d(const type v) { return (v.d_ * type(9.0/20.0)); }
   inline type  notl(const type v) { return (v    != type(0) ? type(0) : type(1)); }
   inline type  frac(const type v) { return (v.d_ - static_cast<long long>(v.d_)); }
   inline type trunc(const type v) { return type((double)static_cast<long long>(v.d_));    }

   inline type modulus(const type v0, const type v1) { return std::fmod(v0.d_,v1.d_);            }
   inline type     pow(const type v0, const type v1) { return std::pow(v0.d_,v1.d_);             }
   inline type    logn(const type v0, const type v1) { return std::log(v0.d_) / std::log(v1.d_); }
   inline type    root(const type v0, const type v1) { return pow(v0,type(1.0) / v1);            }
   inline type   atan2(const type v0, const type v1) { return std::atan2(v0.d_,v1.d_);           }
   inline type     max(const type v0, const type v1) { return std::max(v0.d_,v1.d_);             }
   inline type     min(const type v0, const type v1) { return std::min(v0.d_,v1.d_);             }

   inline bool is_true  (const type v) { return (v != type(0)); }
   inline bool is_false (const type v) { return (v == type(0)); }

   inline type equal(const type v0x, const type v1x)
   {
      const type v0 = v0x.d_;
      const type v1 = v1x.d_;
      static const type epsilon = type(0.0000000001);
      return (abs(v0 - v1) <= (max(type(1),max(abs(v0),abs(v1))) * epsilon)) ? type(1) : type(0);
   }

   inline type expm1(const type vx)
   {
      const type v = vx.d_;
      if (abs(v) < type(0.00001))
         return type(v + (0.5 * v * v));
      else
         return type(exp(v) - type(1));
   }

   inline type nequal(const type v0, const type v1)
   {
      static const type epsilon = type(0.0000000001);
      return (abs(v0 - v1) > (max(type(1),max(abs(v0),abs(v1))) * epsilon)) ? type(1) : type(0);
   }

   inline type log1p(const type v)
   {
      if (v > type(-1))
      {
         if (abs(v) > type(0.0001))
         {
            return log(type(1) + v);
         }
         else
            return (type(-0.5) * v + type(1)) * v;
      }
      else
         return type(std::numeric_limits<double>::quiet_NaN());
   }

   inline type round(const type v)
   {
      return ((v < type(0)) ? ceil(v - type(0.5)) : floor(v + type(0.5)));
   }

   inline type roundn(const type v0, const type v1)
   {
      const type p10 = pow(type(10),trunc(v1));
      if (v0 < type(0))
         return type(ceil ((v0 * p10) - type(0.5)) / p10);
      else
         return type(floor((v0 * p10) + type(0.5)) / p10);
   }

   inline type hypot(const type v0, const type v1)
   {
      return sqrt((v0 * v0) + (v1 * v1));
   }

   inline type shr(const type v0, const type v1)
   {
      return v0 * (type(1) / pow(type(2),trunc(v1)));
   }

   inline type shl(const type v0, const type v1)
   {
      return v0 * pow(type(2),trunc(v1));
   }

   inline type sgn(const type v)
   {
           if (v > type(0)) return type(+1);
      else if (v < type(0)) return type(-1);
      else               return type( 0);
   }

   inline type nand(const type v0, const type& v1)
   {
      return (is_false(v0) || is_false(v1)) ? type(1) : type(0);
   }

   inline type nor(const type v0, const type& v1)
   {
      return (is_false(v0) && is_false(v1)) ? type(1) : type(0);
   }

   inline type xnor(const type v0, const type& v1)
   {
      const bool v0_true = is_true(v0);
      const bool v1_true = is_true(v1);
      if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
         return type(1);
      else
         return type(0);
   }

   inline type erf(type v)
   {
      const type t = type(1) / (type(1) + type(0.5) * abs(v));

      static const type c[] = {
                                type( 1.26551223), type(1.00002368),
                                type( 0.37409196), type(0.09678418),
                                type(-0.18628806), type(0.27886807),
                                type(-1.13520398), type(1.48851587),
                                type(-0.82215223), type(0.17087277)
                              };

      type result = type(1) - t * exp((-v * v) -
                             c[0] + t * (c[1] + t *
                            (c[2] + t * (c[3] + t *
                            (c[4] + t * (c[5] + t *
                            (c[6] + t * (c[7] + t *
                            (c[8] + t * (c[9]))))))))));

      return (v >= type(0)) ? result : -result;
   }

   inline type erfc(type v)
   {
      return type(1) - erf(v);
   }
}

#endif
