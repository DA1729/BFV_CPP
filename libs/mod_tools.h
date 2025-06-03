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

int64_t int_reverse(int64_t a, int n);

template <typename T>
std::vector<T> index_reverse(const std::vector<T>& a, int r);

// polynomial functions

// reduced polynomial multiplication
std::vector<int64_t> red_pol_mul(const std::vector<int64_t>& poly_a, const std::vector<int64_t>& poly_b, int64_t m);

// 2nd version without modulus
std::vector<int64_t> red_pol_mul_2(const std::vector<int64_t>& poly_a, const std::vector<int64_t>& poly_b);

// ntt friendly prime generation
int64_t ntt_friendly_prime(int n, int logq, int lambda, std::mt19937& rng);

// primitive root of unity
bool root_of_unity_check(int64_t w, int64_t m, int64_t q);

std::pair<bool, int64_t> find_primitive_root(int64_t m, int64_t q, std::mt19937& rng);

// bfv parameter generation
std::tuple<int64_t, int64_t, int64_t, int64_t, int64_t>bfv_param_gen(int n, int logq, int lambda, std::mt19937& rng);
#endif