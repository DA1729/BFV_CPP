#ifndef MOD_TOOLS_H
#define MOD_TOOLS_H
#include "prime_utils.h"
#include <cstdint>
#include <random>
#include <tuple>
#include <vector>
#include <utility>
#include <stdexcept>

std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b);

int64_t mod_inv(int64_t a, int64_t m);

int64_t gcd(int64_t a, int64_t b);

uint64_t int_reverse(uint64_t a, int n);

template <typename T>
std::vector<T> index_reverse(const std::vector<T>& a, int r);

// polynomial functions

// reduced polynomial multiplication
std::vector<uint64_t> red_pol_mul(const std::vector<uint64_t>& poly_a, const std::vector<uint64_t>& poly_b, uint64_t m);

// 2nd version without modulus
std::vector<uint64_t> red_pol_mul_2(const std::vector<uint64_t>& poly_a, const std::vector<uint64_t>& poly_b);

// ntt friendly prime generation
uint64_t ntt_friendly_prime(int n, int logq, int lambda, std::mt19937& rng);

// primitive root of unity
bool root_of_unity_check(uint64_t w, uint64_t m, uint64_t q);

std::pair<bool, uint64_t> find_primitive_root(uint64_t m, uint64_t q, std::mt19937& rng);

// bfv parameter generation
std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>bfv_param_gen(int n, int logq, int lambda, std::mt19937& rng);
#endif