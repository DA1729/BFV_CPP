#ifndef PRIME_UTILS_H
#define PRIME_UTILS_H

#include <vector>
#include <cstdint>
#include<random>

int64_t mod_pow(int64_t base, int64_t exp, int64_t mod);
bool miller_rabin(int64_t p, int lambda, std::mt19937& rng);
bool is_prime(int64_t n, int lambda, std::mt19937& rng);
int64_t generate_large_prime(int bit_length, int lambda, std::mt19937& rng);


#endif