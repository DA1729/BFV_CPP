#ifndef POLY_TOOLS_H
#define POLY_TOOLS_H
#include <vector>
#include <cstdint>
#include <random>
#include "mod_tools.h"
#include "prime_utils.h"
#include "ntt_utils.h"


struct ntt_params
{
    std::vector<int64_t> w;       // NTT parameter w
    std::vector<int64_t> w_inv;   // NTT parameter w inverse
    std::vector<int64_t> psi;     // NTT parameter psi
    std::vector<int64_t> psi_inv; // NTT parameter psi inverse

};


struct poly_utils_1
{
    int64_t n;                      // Degree of the polynomial
    int64_t q;                      // Modulus
    std::vector<int64_t> F;         // Coefficients
    ntt_params n_p;                 // ntt parameters
    bool in_ntt;                    // Flag to indicate if in NTT domain

};

poly_utils_1 poly_init(int64_t n, int64_t q, const ntt_params& n_p);

void poly_randomize(poly_utils_1& poly, int64_t B, bool domain = false, int type = 0, double mu = 0.0, double sigma = 1.0);

std::string poly_to_string(const poly_utils_1& poly);

poly_utils_1 poly_add(const poly_utils_1& a, const poly_utils_1& b);

poly_utils_1 poly_sub(const poly_utils_1& a, const poly_utils_1& b);

poly_utils_1 poly_mul(const poly_utils_1& a, const poly_utils_1& b);

poly_utils_1 poly_mod(const poly_utils_1& a, int64_t base);

poly_utils_1 poly_round(const poly_utils_1& a);

bool poly_eq(const poly_utils_1& a, const poly_utils_1& b);

poly_utils_1 poly_neg(const poly_utils_1& a);

poly_utils_1 poly_copy(poly_utils_1& dest, const poly_utils_1& a);

poly_utils_1 to_ntt(const poly_utils_1& a);

poly_utils_1 to_pol(const poly_utils_1& a); 
#endif