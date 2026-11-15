#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <omp.h>

// Hàm kiểm tra một số có phải là số nguyên tố không
inline bool is_prime(long long num) {
    if (num <= 1) return false;
    if (num == 2 || num == 3) return true;
    if (num % 2 == 0 || num % 3 == 0) return false;
    for (long long i = 5; i * i <= num; i += 6) {
        if (num % i == 0 || num % (i + 2) == 0) return false;
    }
    return true;
}

// Đếm số nguyên tố tuần tự trong khoảng [1, n]
long long count_primes_serial(long long n) {
    long long count = 0;
    for (long long i = 1; i <= n; ++i) {
        if (is_prime(i)) {
            count++;
        }
    }
    return count;
}

// Đếm số nguyên tố song song bằng OpenMP với số luồng chỉ định
long long count_primes_parallel(long long n, int num_threads) {
    long long count = 0;
    #pragma omp parallel for reduction(+:count) schedule(guided, 1000) num_threads(num_threads)
    for (long long i = 1; i <= n; ++i) {
        if (is_prime(i)) {
            count++;
        }
    }
    return count;
}

int main(int argc, char* argv[]) {
    long long n = 10000000; // Mặc định 10 triệu
    int num_threads = 4;

    if (argc >= 2) n = std::atoll(argv[1]);
    if (argc >= 3) num_threads = std::atoi(argv[2]);

    std::cout << "Kiem tra dem so nguyen to toi N = " << n << " voi " << num_threads << " luong...\n";

    double t0 = omp_get_wtime();
    long long count = count_primes_parallel(n, num_threads);
    double t1 = omp_get_wtime();

    std::cout << "Ket qua: " << count << " so nguyen to.\n";
    std::cout << "Thoi gian: " << std::fixed << std::setprecision(4) << (t1 - t0) << " giay.\n";

    return 0;
}
