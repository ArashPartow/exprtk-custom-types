#
# **************************************************************
# *         C++ Mathematical Expression Toolkit Library        *
# *                                                            *
# * ExprTk MPFR                                                *
# * Author: Arash Partow                                       *
# * URL: http://www.partow.net/programming/exprtk/index.html   *
# *                                                            *
# * Copyright notice:                                          *
# * Free use of the Mathematical Expression Toolkit Library is *
# * permitted under the guidelines and in accordance with the  *
# * most current version of the MIT License.                   *
# * http://www.opensource.org/licenses/MIT                     *
# *                                                            *
# **************************************************************
#


COMPILER         = -c++
#COMPILER        = -clang++
OPTIMIZATION_OPT = -O1
BASE_OPTIONS     = -pedantic-errors -Wall -Wextra -Werror -Wno-long-long -std=c++11
OPTIONS          = $(BASE_OPTIONS) $(OPTIMIZATION_OPT)
LINKER_OPT       = -L/usr/lib -lstdc++ -lm -lmpfr

BUILD_LIST+=exprtk_mpfr_test
BUILD_LIST+=exprtk_mpfr_benchmark
BUILD_LIST+=exprtk_mpfr_repl

all: $(BUILD_LIST)

$(BUILD_LIST) : %: %.cpp exprtk.hpp
	$(COMPILER) $(OPTIONS) -o $@ $@.cpp $(LINKER_OPT)

pgo: exprtk_mpfr_benchmark.cpp exprtk_mpfr_adaptor.hpp exprtk.hpp
	$(COMPILER) $(BASE_OPTIONS) -O3 -march=native -fprofile-generate -o exprtk_mpfr_benchmark exprtk_mpfr_benchmark.cpp $(LINKER_OPT)
	./exprtk_mpfr_benchmark
	$(COMPILER) $(BASE_OPTIONS) -O3 -march=native -fprofile-use -o exprtk_mpfr_benchmark exprtk_mpfr_benchmark.cpp $(LINKER_OPT)

strip_bin :
	@for f in $(BUILD_LIST); do if [ -f $$f ]; then strip -s $$f; echo $$f; fi done;

valgrind_check:
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes --log-file=exprtk_mpfr_test_valgrind.log -v ./exprtk_mpfr_test
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes --log-file=exprtk_mpfr_benchmark_valgrind.log -v ./exprtk_mpfr_benchmark

clean:
	rm -f core.* *~ *.o *.bak *stackdump gmon.out *.gcda *.gcno *.gcnor *.gch
