#ifndef BFV_LIB_H
#define BFV_LIB_H
#include "./../libs/poly_utils.h"
#include <vector>
#include <cmath>
#include <iostream>

class bfv{
    public:
        int64_t n, q, t, T, l, p;
        double mu, sigma;
        ntt_params qnp;

        poly_utils_1 s_k;
        std::vector<poly_utils_1> p_k;
        std::vector<std::vector<poly_utils_1>> rl_k1;
        std::vector<poly_utils_1> rl_k2;

        bfv(int64_t n, int64_t q, int64_t t, double mu, double sigma, ntt_params qnp);

        std::string to_string() const;
        // main functions

        void secret_key_gen();
        void public_key_gen();
        void eval_key_gen_1();
        void eval_key_gen_2();

        std::vector<poly_utils_1> encryption(const poly_utils_1& m);
        poly_utils_1 decryption(const std::vector<poly_utils_1>& c_t);
        poly_utils_1 decryption_2(const std::vector<poly_utils_1>& c_t);

        std::vector<poly_utils_1> relinearization_1(const std::vector<poly_utils_1>& c_t);
        std::vector<poly_utils_1> relinearization_2(const std::vector<poly_utils_1>& c_t);

        poly_utils_1 int_encode(int64_t m);
        int64_t int_decode(const poly_utils_1& m);

        // homomorphic operations
        std::vector<poly_utils_1> homomorphic_addition(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1);
        std::vector<poly_utils_1> homomorphic_subtraction(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1);
        std::vector<poly_utils_1> homomorphic_multiplication(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1);

        void print_params();

}
;
#endif