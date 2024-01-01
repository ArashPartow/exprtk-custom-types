/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Custom Real type Adaptor                                   *
 * Authors: Arash Partow (1999-2024)                          *
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


#ifndef EXPRTK_REAL_ADAPTOR_HPP
#define EXPRTK_REAL_ADAPTOR_HPP


#include <string>
#include "real_type.hpp"


namespace exprtk
{
   namespace details
   {
      namespace numeric { namespace details
      {
         struct my_real_type_tag;

         template <typename T> inline T const_pi_impl(my_real_type_tag);
         template <typename T> inline T const_e_impl (my_real_type_tag);
      }}

      inline bool is_true (const real::type v);
      inline bool is_false(const real::type v);

      template <typename Iterator>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, real::type& t, details::numeric::details::my_real_type_tag);
   }

   namespace rtl { namespace io
   {
      namespace details
      {
         void print_type(const std::string& fmt, const real::type& v, exprtk::details::numeric::details::my_real_type_tag);
      }
   }}

   using details::is_true;
}

#include "exprtk.hpp"

namespace exprtk
{
   namespace details
   {
      namespace numeric
      {
         namespace details
         {
            struct my_real_type_tag {};

            template<> struct number_type<real::type> { typedef my_real_type_tag type; };

            template <>
            struct epsilon_type<real::type>
            {
               static inline real::type value()
               {
                  const real::type epsilon = real::type(0.000000000100);

                  return epsilon;
               }
            };

            inline bool is_nan_impl(const real::type& v, my_real_type_tag)
            {
               return v.d_ != v.d_;
            }

            template <typename T> inline T   abs_impl(const T v, my_real_type_tag) { return real::abs  (v); }
            template <typename T> inline T  acos_impl(const T v, my_real_type_tag) { return real::acos (v); }
            template <typename T> inline T acosh_impl(const T v, my_real_type_tag) { return real::acosh(v); }
            template <typename T> inline T  asin_impl(const T v, my_real_type_tag) { return real::asin (v); }
            template <typename T> inline T asinh_impl(const T v, my_real_type_tag) { return real::asinh(v); }
            template <typename T> inline T  atan_impl(const T v, my_real_type_tag) { return real::atan (v); }
            template <typename T> inline T atanh_impl(const T v, my_real_type_tag) { return real::atanh(v); }
            template <typename T> inline T  ceil_impl(const T v, my_real_type_tag) { return real::ceil (v); }
            template <typename T> inline T   cos_impl(const T v, my_real_type_tag) { return real::cos  (v); }
            template <typename T> inline T  cosh_impl(const T v, my_real_type_tag) { return real::cosh (v); }
            template <typename T> inline T   exp_impl(const T v, my_real_type_tag) { return real::exp  (v); }
            template <typename T> inline T floor_impl(const T v, my_real_type_tag) { return real::floor(v); }
            template <typename T> inline T   log_impl(const T v, my_real_type_tag) { return real::log  (v); }
            template <typename T> inline T log10_impl(const T v, my_real_type_tag) { return real::log10(v); }
            template <typename T> inline T  log2_impl(const T v, my_real_type_tag) { return real::log2 (v); }
            template <typename T> inline T   neg_impl(const T v, my_real_type_tag) { return -v;             }
            template <typename T> inline T   pos_impl(const T v, my_real_type_tag) { return  v;             }
            template <typename T> inline T   sin_impl(const T v, my_real_type_tag) { return real::sin  (v); }
            template <typename T> inline T  sinh_impl(const T v, my_real_type_tag) { return real::sinh (v); }
            template <typename T> inline T  sqrt_impl(const T v, my_real_type_tag) { return real::sqrt (v); }
            template <typename T> inline T   tan_impl(const T v, my_real_type_tag) { return real::tan  (v); }
            template <typename T> inline T  tanh_impl(const T v, my_real_type_tag) { return real::tanh (v); }
            template <typename T> inline T   cot_impl(const T v, my_real_type_tag) { return real::cot  (v); }
            template <typename T> inline T   sec_impl(const T v, my_real_type_tag) { return real::sec  (v); }
            template <typename T> inline T   csc_impl(const T v, my_real_type_tag) { return real::csc  (v); }
            template <typename T> inline T   r2d_impl(const T v, my_real_type_tag) { return real::r2d  (v); }
            template <typename T> inline T   d2r_impl(const T v, my_real_type_tag) { return real::d2r  (v); }
            template <typename T> inline T   d2g_impl(const T v, my_real_type_tag) { return real::d2g  (v); }
            template <typename T> inline T   g2d_impl(const T v, my_real_type_tag) { return real::g2d  (v); }
            template <typename T> inline T  notl_impl(const T v, my_real_type_tag) { return real::notl (v); }
            template <typename T> inline T  frac_impl(const T v, my_real_type_tag) { return real::frac (v); }
            template <typename T> inline T trunc_impl(const T v, my_real_type_tag) { return real::trunc(v); }

            template <typename T> inline T const_pi_impl(my_real_type_tag) { return real::details::constant::pi; }
            template <typename T> inline T const_e_impl (my_real_type_tag) { return real::details::constant::e;  }

            template <typename T>
            inline int to_int32_impl(const T v, my_real_type_tag)
            {
               return static_cast<int>(v);
            }

            template <typename T>
            inline long long to_int64_impl(const T v, my_real_type_tag)
            {
               return static_cast<long long int>(v);
            }

            template <typename T>
            inline unsigned long long to_uint64_impl(const T v, my_real_type_tag)
            {
               return static_cast<unsigned long long int>(v);
            }

            inline bool is_true_impl (const real::type v)
            {
               return (0.0 != v);
            }

            inline bool is_false_impl(const real::type v)
            {
               return (0.0 == v);
            }

            template <typename T>
            inline T expm1_impl(const T v, my_real_type_tag)
            {
               if (abs(v) < T(0.00001))
                  return v + (T(0.5) * v * v);
               else
                  return exp(v) - T(1);
            }

            template <typename T>
            inline T nequal_impl(const T v0, const T v1, my_real_type_tag)
            {
               const T epsilon  = epsilon_type<T>::value();
               const T eps_norm = (std::max(T(1),std::max(abs_impl(v0,my_real_type_tag()),abs_impl(v1,my_real_type_tag()))) * epsilon);
               return (abs_impl(v0 - v1,my_real_type_tag()) > eps_norm) ? T(1) : T(0);
            }

            template <typename T>
            inline T sgn_impl(const T v, my_real_type_tag)
            {
               if      (v > T(0)) return T(+1);
               else if (v < T(0)) return T(-1);
               else               return T( 0);
            }

            template <typename T>
            inline T log1p_impl(const T v, my_real_type_tag)
            {
               if (v > T(-1))
               {
                  if (abs_impl(v,my_real_type_tag()) > T(0.0001))
                  {
                     return log_impl(T(1) + v,my_real_type_tag());
                  }
                  else
                     return (T(-0.5) * v + T(1)) * v;
               }
               else
                  return T(std::numeric_limits<T>::quiet_NaN());
            }

            template <typename T>
            inline T erf_impl(T v, my_real_type_tag)
            {
               return real::erf(v);
            }

            template <typename T>
            inline T erfc_impl(T v, my_real_type_tag)
            {
               return T(1) - erf_impl(v,my_real_type_tag());
            }

            template <typename T>
            inline T ncdf_impl(T v, my_real_type_tag)
            {
               return T(0.5) * erfc_impl(-(v / T(numeric::constant::sqrt2)),my_real_type_tag());
            }

            template <typename T>
            inline T modulus_impl(const T v0, const T v1, my_real_type_tag)
            {
               return modulus(v0,v1);
            }

            template <typename T>
            inline T pow_impl(const T v0, const T v1, my_real_type_tag)
            {
               return real::pow(v0,v1);
            }

            template <typename T>
            inline T logn_impl(const T v0, const T v1, my_real_type_tag)
            {
               return log(v0) / log(v1);
            }

            template <typename T>
            inline T sinc_impl(T v, my_real_type_tag)
            {
               if (abs_impl(v,my_real_type_tag()) >= std::numeric_limits<T>::epsilon())
                   return(sin_impl(v,my_real_type_tag()) / v);
               else
                  return T(1);
            }

            template <typename T>
            inline T xor_impl(const T v0, const T v1, my_real_type_tag)
            {
               return (is_false_impl(v0) != is_false_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T xnor_impl(const T v0, const T v1, my_real_type_tag)
            {
               const bool v0_true = is_true_impl(v0);
               const bool v1_true = is_true_impl(v1);
               if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
                  return T(1);
               else
                  return T(0);
            }

            template <typename T>
            inline T equal_impl(const T v0, const T v1, my_real_type_tag)
            {
               const T epsilon = epsilon_type<T>::value();
               const T eps_norm = (max(T(1),max(abs_impl(v0,my_real_type_tag()),abs_impl(v1,my_real_type_tag()))) * epsilon);
               return (abs_impl(v0 - v1,my_real_type_tag()) <= eps_norm) ? T(1) : T(0);
            }

            template <typename T>
            inline T round_impl(const T v, my_real_type_tag)
            {
               return ((v < T(0)) ? ceil(v - T(0.5)) : floor(v + T(0.5)));
            }

            template <typename T>
            inline T roundn_impl(const T v0, const T v1, my_real_type_tag)
            {
               const int index = std::max<int>(0, std::min<int>(pow10_size - 1, static_cast<int>(floor(v1))));
               const T p10 = T(pow10[index]);
               if (v0 < T(0))
                  return T(ceil ((v0 * p10) - T(0.5)) / p10);
               else
                  return T(floor((v0 * p10) + T(0.5)) / p10);
            }

            template <typename T>
            inline bool is_integer_impl(const T v, my_real_type_tag)
            {
               return (T(0) == modulus_impl(v,T(1),my_real_type_tag()));
            }

            template <typename T>
            inline T root_impl(const T v0, const T v1, my_real_type_tag)
            {
               return pow(v0,T(1) / v1);
            }

            template <typename T>
            inline T hypot_impl(const T v0, const T v1, my_real_type_tag)
            {
               return sqrt((v0 * v0) + (v1 * v1));
            }

            template <typename T>
            inline T atan2_impl(const T v0, const T v1, my_real_type_tag)
            {
               return std::atan2(v0.d_,v1.d_);
            }

            template <typename T>
            inline T shr_impl(const T v0, const T v1, my_real_type_tag)
            {
               return v0 * (T(1) / pow(T(2),v1));
            }

            template <typename T>
            inline T shl_impl(const T v0, const T v1, my_real_type_tag)
            {
               return v0 * pow(T(2),v1);
            }

            template <typename T>
            inline T and_impl(const T v0, const T v1, my_real_type_tag)
            {
               return (is_true_impl(v0) && is_true_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T nand_impl(const T v0, const T v1, my_real_type_tag)
            {
               return (is_false_impl(v0) || is_false_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T or_impl(const T v0, const T v1, my_real_type_tag)
            {
               return (is_true_impl(v0) || is_true_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T nor_impl(const T v0, const T v1, my_real_type_tag)
            {
               return (is_false_impl(v0) && is_false_impl(v1)) ? T(1) : T(0);
            }
         }
      }

      template <typename Iterator>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, real::type& t, details::numeric::details::my_real_type_tag)
      {
         typename numeric::details::number_type<double>::type num_type;
         return string_to_real<Iterator,double>(itr_external,end,t.d_,num_type);
      }

      inline bool is_true (const real::type v) { return real::is_true(v);  }
      inline bool is_false(const real::type v) { return real::is_false(v); }
   }

   namespace rtl { namespace io
   {
      namespace details
      {
         inline void print_type(const std::string& fmt, const real::type& v, exprtk::details::numeric::details::my_real_type_tag)
         {
            #if defined(__clang__)
               #pragma clang diagnostic push
               #pragma clang diagnostic ignored "-Wformat-nonliteral"
            #elif defined(__GNUC__) || defined(__GNUG__)
               #pragma GCC diagnostic push
               #pragma GCC diagnostic ignored "-Wformat-nonliteral"
            #elif defined(_MSC_VER)
            #endif

            printf(fmt.c_str(),static_cast<double>(v));

            #if defined(__clang__)
               #pragma clang diagnostic pop
            #elif defined(__GNUC__) || defined(__GNUG__)
               #pragma GCC diagnostic pop
            #elif defined(_MSC_VER)
            #endif
         }
      }
   }}
}

#endif
