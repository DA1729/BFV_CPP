#ifndef PRIME_UTILS_H
#define PRIME_UTILS_H

#include <vector>
#include <cstdint>
#include<random>

uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod);
bool miller_rabin(uint64_t p, int lambda, std::mt19937& rng);
bool is_prime(uint64_t n, int lambda, std::mt19937& rng);
uint64_t generate_large_prime(int bit_length, int lambda, std::mt19937& rng);


#endif