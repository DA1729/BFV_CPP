#include "prime_utils.h"
#include <cmath>
#include <stdexcept>
#include <ctime>

uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1)
            result = (__uint128_t(result) * base) % mod;
        base = (__uint128_t(base) * base) % mod;
        exp >>= 1;
    }
    return result;
}

std::vector<int> get_low_primes() {
    return {
        3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
        101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
        191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
        281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
        389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
        491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
        607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
        719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
        829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
        953,967,971,977,983,991,997
    };
}

bool miller_rabin(uint64_t p, int lambda, std::mt19937& rng) {
    uint64_t r = p - 1;
    int u = 0;
    while ((r & 1) == 0) {
        ++u;
        r >>= 1;
    }

    std::uniform_int_distribution<uint64_t> dist(2, p - 2);
    for (int i = 0; i < lambda; ++i) {
        uint64_t a = dist(rng);
        uint64_t z = mod_pow(a, r, p);
        if (z != 1 && z != p - 1) {
            for (int j = 0; j < u - 1; ++j) {
                z = mod_pow(z, 2, p);
                if (z == 1) return false;
                if (z == p - 1) break;
            }
            if (z != p - 1) return false;
        }
    }

    return true;
}

bool is_prime(uint64_t n, int lambda, std::mt19937& rng) {
    if (n < 2) return false;
    if (n == 2) return true;

    auto low_primes = get_low_primes();
    for (int p : low_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    return miller_rabin(n, lambda, rng);
}

uint64_t generate_large_prime(int bit_length, int lambda, std::mt19937& rng) {
    int attempts = 100 * (std::log2(bit_length) + 1);
    std::uniform_int_distribution<uint64_t> dist((1ULL << (bit_length - 1)), (1ULL << bit_length) - 1);

    while (attempts-- > 0) {
        uint64_t candidate = dist(rng);
        if (is_prime(candidate, lambda, rng)) {
            return candidate;
        }
    }

    throw std::runtime_error("Failed to generate a large prime after many attempts");
}
