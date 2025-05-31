#ifndef NTT_UTILS_H
#define NTT_UTILS_H
#include "mod_tools.h"
#include <cstdint>
#include <random>
#include <tuple>
#include <vector>
#include <utility>
#include <stdexcept>

std::vector<uint64_t> ntt(const std::vector<uint64_t>& input, const std::vector<uint64_t>& w_table, uint64_t q);

std::vector<uint64_t> intt(const std::vector<uint64_t>& input, const std::vector<uint64_t>& w_table, uint64_t q);
#endif