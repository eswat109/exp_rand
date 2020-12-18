#include <iostream>
#include <omp.h>
#include <mutex>
#include <cmath>
#include <fstream>

#include <cstdio>
#include <vector>
#include <thread>
#include <ctime>

#define ElemCount 100000000
#define CACHE_LINE 64
#define A 1103515245u
#define b 12345u
#define c 4294967296u

using namespace std;

void run_experiment_random(double(*f)(std::vector<unsigned>, unsigned, unsigned, unsigned, unsigned)) {

    //ofstream cout("rand.txt");

    unsigned T = omp_get_num_procs();
    double time, res;
    //unsigned* arr = (unsigned*)malloc(ElemCount * sizeof(unsigned));
    std::vector<unsigned> arr;
    arr.resize(ElemCount);
    cout << "Threads:\t" << "\tresult:\t" << "\ttime:\t" << "\tspeedup:" << "\n";
    for (int t = 1; t <= T; t++) {
        omp_set_num_threads(t);
        double time1 = omp_get_wtime();
        res = f(arr, ElemCount, 1, 100, t);
        double time2 = omp_get_wtime();
        if (t == 1) {
            time = time2 - time1;
        }
        cout << t << "\t\t" << res << "\t\t" << (time2 - time1) << "\t\t" << (time / (time2 - time1)) << "\n";
    }
    cout << "\n";

    //cout.close();

}



std::vector<unsigned> get_A(unsigned T) {

    std::vector<unsigned> result;

    result.reserve(T);

    result.emplace_back(A);

    for (unsigned i = 1; i < T + 1; i++) {

        result.emplace_back((result.back() * A) % c);

    }

    return result;

}


struct aligned_unsigned {
    alignas(CACHE_LINE) unsigned value;
};

double randomize(std::vector<unsigned> V, unsigned N, unsigned mina, unsigned maxa, unsigned nt) {

    unsigned T = nt;

    std::vector<unsigned> multipliers = get_A(T);

    double median = 0.0;

    std::mutex mtx;

    std::vector<aligned_unsigned> partial_rand(T);

    std::vector<std::thread> threads;

    unsigned seed = std::time(nullptr);

    for (std::size_t t = 0; t < T; ++t) {
        threads.emplace_back([t, T, &V, N, seed, &multipliers, &mtx, &median, &mina, &maxa]() {

            unsigned _A = multipliers.back();
            unsigned off = (b * (_A - 1) / (A - 1)) % c;

            unsigned x = ((seed * multipliers[t]) % c + (b * (multipliers[t] - 1) / (A - 1)) % c) % c;

            double my_median = 0.0;

            for (size_t i = t; i < N; i += T) {

                V[i] = (x % (maxa - mina)) + mina;

                my_median += (double)V[i];

                x = ((x * _A) % c + off) % c;

            }
            mtx.lock();
            median += my_median;
            mtx.unlock();

            });
    }
    for (auto& thread : threads) { thread.join(); }
    return median / N;
}

int main() {

    //run_experiment_random(omp_auto_reduction);

    run_experiment_random(randomize);

    return 0;
}