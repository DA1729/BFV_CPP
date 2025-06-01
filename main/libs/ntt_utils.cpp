#include <stdexcept>
#include <cmath>
#include <ctime>
#include "ntt_utils.h"

std::vector<int64_t> ntt(const std::vector<int64_t>& input, const std::vector<int64_t>& w_table, int64_t q){
    size_t n = input.size();
    std::vector<int64_t> output(input);

    int levels = std::log2(n);
    
    if ((1ULL << levels) != n){
        throw std::invalid_argument("input size must be a power of two");
    }

    for (int i = 0; i < levels; i++){
        int m = 1 << i;
        int step = 1 << (levels - i - 1);
        
        for (int j = 0; j < m; j++){
            for (int k = 0; k < step; k++){
                int s = j*(step << 1) + k;
                int t = s + step;

                int64_t w = w_table[(1 << i) *k];
                int64_t u = output[s];
                int64_t v = output[t];

                output[s] = (u + v) % q;
                output[t] = ((q + u - v)*w)%q;
            }
        }
    }

    output = index_reverse(output, levels);
    return output;
}

std::vector<int64_t> intt(const std::vector<int64_t>& input, const std::vector<int64_t>& w_table, int64_t q){
    size_t n = input.size();
    std::vector<int64_t> output(input);

    int levels = std::log2(n);
    if ((1ULL << levels) != n) {
        throw std::invalid_argument("Input size must be a power of two.");
    }

    for (int i = 0; i < levels; ++i) {
        int m = 1 << i;
        int step = 1 << (levels - i - 1);
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < step; ++k) {
                int s = j * (step << 1) + k;
                int t = s + step;

                int64_t w = w_table[(1 << i) * k];
                int64_t u = output[s];
                int64_t v = output[t];

                output[s] = (u + v) % q;
                output[t] = ((q + u - v) * w) % q;
            }
        }
    }

    output = index_reverse(output, levels);

    int64_t n_inv = mod_inv(n, q);
    for (auto& val : output) {
        val = (val * n_inv) % q;
    }

    return output;
}