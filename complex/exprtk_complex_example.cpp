/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Example using a simple Complex type                        *
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


#include <cstdio>
#include <string>
#include "exprtk_complex_adaptor.hpp"
#include "exprtk.hpp"


template <typename T>
void complex_numbers()
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>   expression_t;
   typedef exprtk::parser<T>       parser_t;

   const std::vector<std::string> expressions =
   {
      "(1 + i) / (3 + 2i)",
      "x + y",
      "x / y",
      "(x + i) / (3y + 2i)",
   };

   symbol_table_t symbol_table;

   T i = T(0.0,1.0);
   T x = T(1.1,0.0);
   T y = T(2.2,0.0);
   T v[5];

   symbol_table.add_variable ("i" , i );
   symbol_table.add_variable ("x" , x );
   symbol_table.add_variable ("y" , y );
   symbol_table.add_vector   ("v" , v );

   for (std::size_t i = 0; i < expressions.size(); ++i)
   {
      expression_t expression;
      expression.register_symbol_table(symbol_table);

      parser_t parser;
      if (!parser.compile(expressions[i], expression))
      {
         for (std::size_t error_idx = 0; error_idx < parser.error_count(); ++error_idx)
         {
            const auto error = parser.get_error(error_idx);
            printf("Err: %02d Pos: %02d Type: [%14s] Msg: %s\tExpression: %s\n",
                   static_cast<unsigned int>(error_idx),
                   static_cast<unsigned int>(error.token.position),
                   exprtk::parser_error::to_str(error.mode).c_str(),
                   error.diagnostic.c_str(),
                   expressions[i].c_str());
         }
      }

      const auto result = expression.value();

      printf("%s = %12.5f + %12.5fi\n",
             expressions[i].c_str(),
             result.c_.real(), result.c_.imag());
   }

   return;
}

int main()
{
   complex_numbers<cmplx::complex_t>();
   return 0;
}
