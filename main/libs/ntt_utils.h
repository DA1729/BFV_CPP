#ifndef NTT_UTILS_H
#define NTT_UTILS_H
#include "mod_tools.h"
#include <cstdint>
#include <random>
#include <tuple>
#include <vector>
#include <utility>
#include <stdexcept>

std::vector<int64_t> ntt(const std::vector<int64_t>& input, const std::vector<int64_t>& w_table, int64_t q);

std::vector<int64_t> intt(const std::vector<int64_t>& input, const std::vector<int64_t>& w_table, int64_t q);

template <typename T>
std::vector<T> index_reverse(const std::vector<T>& a, int r){
    size_t n = a.size();
    std::vector<T> b(n);

    for (size_t i = 0; i < n; i++){
        size_t index = int_reverse(i, r);
        b[index] = a[i];

    }

    return b;
}

#endif