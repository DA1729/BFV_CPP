#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include "./bfv_lib/bfv_lib.h"
#include "./libs/mod_tools.h"

int main() {
    std::cout << "--- Starting BFV Demo" << std::endl;

    // Parameter Setup
    int PD = 0;  // 0: generate -- 1: pre-defined

    uint64_t t, n, q, psi, psiv, w, wv;

    if (PD == 0) {
        t = 16;
        n = 1024;
        q = 132120577;
        psi = 73993;

        psiv = mod_inv(psi, q);
        w    = (psi * psi) % q;
        wv   = mod_inv(w, q);
    } else {
        t = 16; n = 1024;
        int logq = 27;
        std::mt19937 rng(time(0));
        std::tie(q, psi, psiv, w, wv) = bfv_param_gen(n, logq, 128, rng);
    }

    double mu = 0;
    double sigma = 0.5 * 3.2;
    uint64_t T = 256;
    uint64_t p = q*q*q + 1;

    // Table Generation
    std::vector<uint64_t> w_table(n), wv_table(n), psi_table(n), psiv_table(n);
    w_table[0] = wv_table[0] = psi_table[0] = psiv_table[0] = 1;

    for (size_t i = 1; i < n; i++) {
        w_table[i] = (w_table[i - 1] * w) % q;
        wv_table[i] = (wv_table[i - 1] * wv) % q;
        psi_table[i] = (psi_table[i - 1] * psi) % q;
        psiv_table[i] = (psiv_table[i - 1] * psiv) % q;
    }

    std::vector<std::vector<uint64_t>> qnp = {w_table, wv_table, psi_table, psiv_table};

    bfv Evaluator(n, q, t, T, 0, p, mu, sigma, qnp);

    // Key Generation
    Evaluator.secret_key_gen();
    Evaluator.public_key_gen();
    Evaluator.eval_key_gen_1(T);
    Evaluator.eval_key_gen_2(p);

    // Random messages
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int64_t> dist(-(1 << 15), (1 << 15) - 1);

    int64_t n1 = dist(rng);
    int64_t n2 = dist(rng);

    std::cout << "--- Random integers generated:\n* n1: " << n1 << "\n* n2: " << n2;
    std::cout << "\n* n1+n2: " << (n1 + n2) << "\n* n1-n2: " << (n1 - n2);
    std::cout << "\n* n1*n2: " << (n1 * n2) << "\n";

    // Encoding
    auto m1 = Evaluator.int_encode(n1);
    auto m2 = Evaluator.int_encode(n2);

    // Encryption
    auto ct1 = Evaluator.encryption(m1);
    auto ct2 = Evaluator.encryption(m2);

    // Homomorphic Addition
    auto ct_add = Evaluator.homomorphic_addition(ct1, ct2);
    auto dec_add = Evaluator.decryption(ct_add);
    int64_t res_add = Evaluator.int_decode(dec_add);

    std::cout << "\n--- Homomorphic Addition Result: " << res_add << std::endl;
    std::cout << (res_add == (n1 + n2) ? "* Homomorphic addition works.\n" : "* Homomorphic addition FAILED.\n");

    // Homomorphic Subtraction
    auto ct_sub = Evaluator.homomorphic_subtraction(ct1, ct2);
    auto dec_sub = Evaluator.decryption(ct_sub);
    int64_t res_sub = Evaluator.int_decode(dec_sub);

    std::cout << "\n--- Homomorphic Subtraction Result: " << res_sub << std::endl;
    std::cout << (res_sub == (n1 - n2) ? "* Homomorphic subtraction works.\n" : "* Homomorphic subtraction FAILED.\n");

    // Homomorphic Multiplication (no relinearization)
    auto ct_mul = Evaluator.homomorphic_multiplication(ct1, ct2);
    auto dec_mul = Evaluator.decryption_2(ct_mul);
    int64_t res_mul = Evaluator.int_decode(dec_mul);

    std::cout << "\n--- Homomorphic Multiplication Result (no relinearization): " << res_mul << std::endl;
    std::cout << (res_mul == (n1 * n2) ? "* Homomorphic multiplication works.\n" : "* Homomorphic multiplication FAILED.\n");

    // Homomorphic Multiplication + Relinearization v1
    auto ct_mul_rel1 = Evaluator.homomorphic_multiplication(ct1, ct2);
    ct_mul_rel1 = Evaluator.relinearization_1(ct_mul_rel1);
    auto dec_mul_rel1 = Evaluator.decryption(ct_mul_rel1);
    int64_t res_mul_rel1 = Evaluator.int_decode(dec_mul_rel1);

    std::cout << "\n--- Homomorphic Multiplication + Relinearization v1 Result: " << res_mul_rel1 << std::endl;
    std::cout << (res_mul_rel1 == (n1 * n2) ? "* Relinearized multiplication v1 works.\n" : "* Relinearized multiplication v1 FAILED.\n");

    return 0;
}
