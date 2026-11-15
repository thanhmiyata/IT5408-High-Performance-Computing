CXX ?= clang++
MPICXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
THREAD_FLAGS = -pthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BREW_PREFIX := $(shell brew --prefix libomp 2>/dev/null || echo "/opt/homebrew/opt/libomp")
    OMP_FLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/include -L$(BREW_PREFIX)/lib -lomp
else
    OMP_FLAGS = -fopenmp
endif

TARGETS = tsp_solver prime_count count_paths matmul_serial matmul_thread matmul_openmp matmul_mpi

all: $(TARGETS)

tsp_solver: tsp_solver.cpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $< -o $@

prime_count: prime_count.cpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $< -o $@

count_paths: count_paths.cpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $< -o $@

matmul_serial: matmul_serial.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

matmul_thread: matmul_thread.cpp common.hpp
	$(CXX) $(CXXFLAGS) $(THREAD_FLAGS) $< -o $@

matmul_openmp: matmul_openmp.cpp common.hpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $< -o $@

matmul_mpi: matmul_mpi.cpp common.hpp
	$(MPICXX) $(CXXFLAGS) $< -o $@

test: all
	@echo "=== TEST TSP (tsp7.inp) ==="
	./tsp_solver tsp7.inp 4
	@echo "\n=== TEST PRIME COUNT (N=1,000,000) ==="
	./prime_count 1000000 4
	@echo "\n=== TEST COUNT PATHS (graph01.txt, k=3) ==="
	./count_paths graph01.txt 3 4
	@echo "\n=== TEST MATMUL SERIAL (512x512x512) ==="
	./matmul_serial 512 512 512
	@echo "\n=== TEST MATMUL OPENMP (512x512x512, 4T) ==="
	./matmul_openmp 512 512 512 4

benchmark: all
	source .venv/bin/activate && python3 benchmark_week2.py

pdf:
	source .venv/bin/activate && python3 export_pdf.py

zip: pdf
	@echo "File ZIP Bai_Tap_Tuan_2_20242422M.zip da duoc dong goi san sang."

clean:
	rm -f $(TARGETS) *.o *.dSYM Bao_cao_Bai_tap_Tuan_2.html

.PHONY: all test benchmark pdf zip clean
