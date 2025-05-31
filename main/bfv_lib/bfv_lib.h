#ifndef BFV_LIB_H
#define BFV_LIB_H
#include <vector>
#include <cstdint>
#include <cstddef>
#include <random>
#include <chrono>
#include <iostream>
#include <utility>
#include <string>
#include "./../libs/mod_tools.h"
#include "./../libs/ntt_utils.h"
#include "./../libs/prime_utils.h"
#include "./../libs/poly_tools.h"


//using poly_tools = std::vector<int64_t>;

class bfv {
public:
    size_t n;
    uint64_t q, t, T, l, p;
    double mu, sigma;
    std::vector<std::vector<uint64_t>> qnp;

    poly_tools s_k;  
    std::vector<poly_tools> p_k;
    std::vector<std::pair<poly_tools, poly_tools>> rl_k_1;
    std::vector<poly_tools> rl_k_2;

    bfv(size_t n, uint64_t q, uint64_t t, uint64_t T, uint64_t l, uint64_t p,
        double mu, double sigma, const std::vector<std::vector<uint64_t>>& qnp);

    void secret_key_gen();
    void public_key_gen();
    void eval_key_gen_1(uint64_t T_in);
    void eval_key_gen_2(uint64_t p_in);

    std::vector<poly_tools> encryption(const poly_tools& m_in);
    poly_tools decryption(const std::vector<poly_tools> &c_t_in);
    poly_tools decryption_2(const std::vector<poly_tools> &c_t_in);

    // relinearization
    std::vector<poly_tools> relinearization_1(const std::vector<poly_tools> &ct);
    std::vector<poly_tools> relinearization_2(const std::vector<poly_tools> &ct);


    // encoding and decoding
    poly_tools int_encode(int64_t m);
    int64_t int_decode(const poly_tools& m);


    // homomorphic evaluations
    std::vector<poly_tools> homomorphic_addition(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1);
    std::vector<poly_tools> homomorphic_subtraction(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1);
    std::vector<poly_tools> homomorphic_multiplication(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1);

};

#endif
