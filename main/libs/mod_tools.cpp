#include "mod_tools.h"
#include "prime_utils.h"
#include <stdexcept>
#include <cmath>
#include <ctime>


std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b){
    if (a == 0){
        return {b, 0, 1};
    }

    else{
        auto [g, y, x] = extended_gcd(b%a, a);
        return {g, x-(b/a)*y, y};
    }
}
//  modular inverse 
int64_t mod_inv(int64_t a, int64_t m){
    auto[g, x, y] = extended_gcd(a, m);

    if (g != 1){
        throw std::invalid_argument("modular inverse doesn't exist");
    }

    else{
        return (x % m + m)%m;
    }
}


// gcd
int64_t gcd(int64_t a, int64_t b){
    int64_t a_1 = a;
    int64_t b_1 = b;

    while (b_1 != 0){
        int64_t temp = b_1;
        b_1 = a_1%b_1;
        a_1 = temp;
    }

    return a_1;
}

uint64_t int_reverse(uint64_t a, int n){
    uint64_t result = 0;

    for (int i = 0; i < n; i++){
        result <<= 1;
        result |= (a&1);
        a >>= 1;
    }

    return result;
}


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


// polynomial functions

// reduced polynomial multiplication
std::vector<uint64_t> red_pol_mul(const std::vector<uint64_t>& poly_a, const std::vector<uint64_t>& poly_b, uint64_t m){
    size_t degree = poly_a.size();

    std::vector<int64_t> product_coeffs(2*degree, 0);
    std::vector<uint64_t> reduced_poly(degree, 0);

    // normal multiplication modulo m
    for (size_t i = 0; i < degree; i++){
        for (size_t j = 0; j < degree; j++){
            product_coeffs[i + j] = (product_coeffs[i + j] + poly_a[i] * poly_b[j]) % m;

        }
    }

    // reducing modulo x^n + 1
    for (size_t i = 0; i < degree; i++){
        reduced_poly[i] = (product_coeffs[i] - product_coeffs[i + degree])%m;
        if (reduced_poly[i] < 0){
            reduced_poly[i] += m;
        }
    }

    return reduced_poly;
}

std::vector<uint64_t> red_pol_mul_2(const std::vector<uint64_t>& poly_a, const std::vector<uint64_t>& poly_b){
    size_t degree = poly_a.size();

    std::vector<int64_t> product_coeffs(2*degree, 0);
    std::vector<uint64_t> reduced_poly(degree, 0);

    // normal multiplication modulo m
    for (size_t i = 0; i < degree; i++){
        for (size_t j = 0; j < degree; j++){
            product_coeffs[i + j] = (product_coeffs[i + j] + poly_a[i] * poly_b[j]);

        }
    }

    // reducing modulo x^n + 1
    for (size_t i = 0; i < degree; i++){
        reduced_poly[i] = (product_coeffs[i] - product_coeffs[i + degree]);
    }

    return reduced_poly;
}


// ntt friendly prime generation
uint64_t ntt_friendly_prime(int n, int logq, int lambda, std::mt19937& rng){
    // the log is base 2
    int64_t step = 2 * n;
    int64_t candidate = (1ULL << logq) - step + 1;
    int64_t l_bound = 1ULL << (logq - 1);

    while (candidate > l_bound){
        if (is_prime(candidate, lambda, rng)){
            return candidate;
        }
        else{
            candidate -= step;
        }
    }

    throw std::runtime_error("failed to find a suitable ntt prime");
}


// primite root of unity
bool root_of_unity_check(uint64_t w, uint64_t m, uint64_t q){
    if (w == 0){
        return false;
    }

    else {
        uint64_t power = mod_pow(w, m/2, q);
        if (power == (q-1)){
            return true;
        }

        else{
            return false;
        }
    }
}

std::pair<bool, uint64_t> find_primitive_root(uint64_t m, uint64_t q, std::mt19937& rng){
    uint64_t g = (q - 1)/m;

    if ((q-1) != g*m){                  // cheacking if m divides q-1 exactly
        return std::make_pair(false, 0);
    }

    std::uniform_int_distribution<uint64_t> dist(2, q-1);

    int count = 0;
    const int max_count = 100;

    while (count < max_count){
        uint64_t a = dist(rng);

        uint64_t b = mod_pow(a, g, q);

        if (root_of_unity_check(b, m, q)){
            return std::make_pair(true, b);
        }
        count++;
    }

    return std::make_pair(true, 0);

}

// bfv parameter generation 
std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> bfv_param_gen(int n, int logq, int lambda, std::mt19937& rng){
    bool p_found = false;
    uint64_t q;
    uint64_t psi;

    while (!p_found){
        q = ntt_friendly_prime(n, logq, lambda, rng);
        std::tie(p_found, psi) = find_primitive_root(2*n, q, rng);

    }

    uint64_t psi_inv = mod_inv(psi, q);
    uint64_t w = mod_pow(psi, 2, q);
    uint64_t w_inv = mod_inv(w, q);

    return std::make_tuple(q, psi, psi_inv, w, w_inv);
}