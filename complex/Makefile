#
# **************************************************************
# *         C++ Mathematical Expression Toolkit Library        *
# *                                                            *
# * ExprTk Simple Complex type                                 *
# * Author: Arash Partow (1999-2025)                           *
# * URL: https://www.partow.net/programming/exprtk/index.html  *
# *                                                            *
# * Copyright notice:                                          *
# * Free use of the Mathematical Expression Toolkit Library is *
# * permitted under the guidelines and in accordance with the  *
# * most current version of the MIT License.                   *
# * https://www.opensource.org/licenses/MIT                    *
# * SPDX-License-Identifier: MIT                               *
# *                                                            *
# **************************************************************
#


COMPILER         = -c++
#COMPILER        = -clang++
OPTIMIZATION_OPT = -O1 -DNDEBUG
BASE_OPTIONS     = -pedantic-errors -Wall -Wextra -Werror -std=c++17
OPTIONS          = $(BASE_OPTIONS) $(OPTIMIZATION_OPT)
LINKER_OPT       = -L/usr/lib -lstdc++ -lm

BUILD_LIST+=exprtk_complex_example

all: $(BUILD_LIST)

$(BUILD_LIST) : %: %.cpp exprtk.hpp
	$(COMPILER) $(OPTIONS) -o $@ $@.cpp $(LINKER_OPT)

clean:
	rm -f core.* *~ *.o *.bak *stackdump gmon.out *.gcda *.gcno *.gcnor *.gch
