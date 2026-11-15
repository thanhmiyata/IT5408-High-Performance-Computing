CXX = clang++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

# OpenMP flags cho macOS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BREW_PREFIX := $(shell brew --prefix libomp 2>/dev/null || echo "/opt/homebrew/opt/libomp")
    OMP_FLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/include -L$(BREW_PREFIX)/lib -lomp
else
    OMP_FLAGS = -fopenmp
endif

all: prime_count

prime_count: prime_count.cpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $< -o $@

run: prime_count
	./prime_count 10000000 4

benchmark: prime_count
	source .venv/bin/activate && python3 run_benchmark.py

pdf:
	source .venv/bin/activate && python3 export_pdf.py

zip: pdf
	@echo "File ZIP Bai_Tap_Tuan_1_20242422M.zip da duoc tao san sang de nop form."

clean:
	rm -f prime_count *.o *.dSYM Bao_cao_Bai_tap_Tuan_1.html

.PHONY: all run benchmark pdf zip clean
