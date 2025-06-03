#include "./bfv_lib/bfv_lib.h"
#include "./libs/mod_tools.h"
#include <iostream>
#include <random>
#include <cmath>

int main() {
    std::cout << "--- Starting BFV Demo" << std::endl;
    
    // Parameter generation (corresponding to Python code)
    int64_t t = 16;
    int64_t n = 1024;
    int64_t q = 132120577;  // log(q) = 27
    int64_t psi = 73993;
    
    // Calculate other necessary parameters
    int64_t psiv = mod_inv(psi, q);
    int64_t w = mod_pow(psi, 2, q);
    int64_t wv = mod_inv(w, q);
    
    // Determine mu, sigma (for discrete gaussian distribution)
    double mu = 0.0;
    double sigma = 0.5 * 3.2;
    
    // Generate polynomial arithmetic tables
    std::vector<int64_t> w_table(n, 1);
    std::vector<int64_t> wv_table(n, 1);
    std::vector<int64_t> psi_table(n, 1);
    std::vector<int64_t> psiv_table(n, 1);
    
    for (int64_t i = 1; i < n; i++) {
        w_table[i] = (w_table[i-1] * w) % q;
        wv_table[i] = (wv_table[i-1] * wv) % q;
        psi_table[i] = (psi_table[i-1] * psi) % q;
        psiv_table[i] = (psiv_table[i-1] * psiv) % q;
    }
    
    // Create ntt_params structure
    ntt_params qnp;
    qnp.w = w_table;
    qnp.w_inv = wv_table;
    qnp.psi = psi_table;
    qnp.psi_inv = psiv_table;
    
    // Generate BFV evaluator
    bfv Evaluator(n, q, t, mu, sigma, qnp);
    
    // Generate Keys
    Evaluator.secret_key_gen();
    Evaluator.public_key_gen();
    Evaluator.eval_key_gen_1();
    // Note: eval_key_gen_2() requires setting p parameter first
    // Evaluator.p = q*q*q + 1; // Set p for relinearization v2
    // Evaluator.eval_key_gen_2();
    
    // Print system parameters
    Evaluator.print_params();
    
    // Generate random messages
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dis(-(1LL << 15), (1LL << 15) - 1);
    
    int64_t n1 = dis(gen);
    int64_t n2 = dis(gen);
    
    std::cout << "\n--- Random integers n1 and n2 are generated." << std::endl;
    std::cout << "* n1: " << n1 << std::endl;
    std::cout << "* n2: " << n2 << std::endl;
    std::cout << "* n1+n2: " << (n1 + n2) << std::endl;
    std::cout << "* n1-n2: " << (n1 - n2) << std::endl;
    std::cout << "* n1*n2: " << (n1 * n2) << std::endl;
    std::cout << std::endl;
    
    // Encode random messages into plaintext polynomials
    std::cout << "--- n1 and n2 are encoded as polynomials m1(x) and m2(x)." << std::endl;
    poly_utils_1 m1 = Evaluator.int_encode(n1);
    poly_utils_1 m2 = Evaluator.int_encode(n2);
    
    std::cout << "* m1(x): [";
    for (int64_t i = 0; i < 10; i++) { // Print first 10 coefficients
        std::cout << m1.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "* m2(x): [";
    for (int64_t i = 0; i < 10; i++) { // Print first 10 coefficients
        std::cout << m2.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    std::cout << std::endl;
    
    // Encrypt messages
    std::vector<poly_utils_1> ct1 = Evaluator.encryption(m1);
    std::vector<poly_utils_1> ct2 = Evaluator.encryption(m2);
    
    std::cout << "--- m1 and m2 are encrypted as ct1 and ct2." << std::endl;
    std::cout << "* ct1[0]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct1[0].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "* ct1[1]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct1[1].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    std::cout << std::endl;
    
    // Test basic encryption/decryption without homomorphic operations
    std::cout << "--- Testing basic encryption/decryption (no homomorphic operations)" << std::endl;
    
    // Decrypt ct1 and ct2 to verify they match original messages
    poly_utils_1 decrypted_m1 = Evaluator.decryption(ct1);
    poly_utils_1 decrypted_m2 = Evaluator.decryption(ct2);
    
    int64_t decoded_n1 = Evaluator.int_decode(decrypted_m1);
    int64_t decoded_n2 = Evaluator.int_decode(decrypted_m2);
    
    std::cout << "* Original n1: " << n1 << ", Decrypted n1: " << decoded_n1 << std::endl;
    std::cout << "* Original n2: " << n2 << ", Decrypted n2: " << decoded_n2 << std::endl;
    
    if (decoded_n1 == n1 && decoded_n2 == n2) {
        std::cout << "* Basic encryption/decryption works correctly." << std::endl;
    } else {
        std::cout << "* Basic encryption/decryption failed." << std::endl;
    }
    std::cout << std::endl;
    
    // Homomorphic Addition
    std::vector<poly_utils_1> ct_add = Evaluator.homomorphic_addition(ct1, ct2);
    poly_utils_1 mt_add = Evaluator.decryption(ct_add);
    int64_t nr_add = Evaluator.int_decode(mt_add);
    int64_t ne_add = n1 + n2;
    
    std::cout << "--- Performing ct_add = Enc(m1) + Enc(m2)" << std::endl;
    std::cout << "* ct_add[0]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct_add[0].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dec = Dec(ct_add)" << std::endl;
    std::cout << "* ct_dec: [";
    for (int64_t i = 0; i < 10; i++) {
        std::cout << mt_add.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dcd = Decode(ct_dec)" << std::endl;
    std::cout << "* ct_dcd: " << nr_add << std::endl;
    
    if (nr_add == ne_add) {
        std::cout << "* Homomorphic addition works." << std::endl;
    } else {
        std::cout << "* Homomorphic addition does not work." << std::endl;
    }
    std::cout << std::endl;
    
    // Homomorphic Subtraction
    std::vector<poly_utils_1> ct_sub = Evaluator.homomorphic_subtraction(ct1, ct2);
    poly_utils_1 mt_sub = Evaluator.decryption(ct_sub);
    int64_t nr_sub = Evaluator.int_decode(mt_sub);
    int64_t ne_sub = n1 - n2;
    
    std::cout << "--- Performing ct_sub = Enc(m1) - Enc(m2)" << std::endl;
    std::cout << "* ct_sub[0]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct_sub[0].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dec = Dec(ct_sub)" << std::endl;
    std::cout << "* ct_dec: [";
    for (int64_t i = 0; i < 10; i++) {
        std::cout << mt_sub.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dcd = Decode(ct_dec)" << std::endl;
    std::cout << "* ct_dcd: " << nr_sub << std::endl;
    
    if (nr_sub == ne_sub) {
        std::cout << "* Homomorphic subtraction works." << std::endl;
    } else {
        std::cout << "* Homomorphic subtraction does not work." << std::endl;
    }
    std::cout << std::endl;
    
    // Multiply two messages (no relinearization)
    std::vector<poly_utils_1> ct_mul = Evaluator.homomorphic_multiplication(ct1, ct2);
    poly_utils_1 mt_mul = Evaluator.decryption_2(ct_mul); // Use decryption_2 for 3-element ciphertext
    int64_t nr_mul = Evaluator.int_decode(mt_mul);
    int64_t ne_mul = n1 * n2;
    
    std::cout << "--- Performing ct_mul = Enc(m1) * Enc(m2) (no relinearization)" << std::endl;
    std::cout << "* ct_mul[0]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct_mul[0].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dec = Dec(ct_mul)" << std::endl;
    std::cout << "* ct_dec: [";
    for (int64_t i = 0; i < 10; i++) {
        std::cout << mt_mul.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dcd = Decode(ct_dec)" << std::endl;
    std::cout << "* ct_dcd: " << nr_mul << std::endl;
    
    if (nr_mul == ne_mul) {
        std::cout << "* Homomorphic multiplication works." << std::endl;
    } else {
        std::cout << "* Homomorphic multiplication does not work." << std::endl;
    }
    std::cout << std::endl;
    
    // Multiply two messages (relinearization v1)
    std::vector<poly_utils_1> ct_mul_relin = Evaluator.homomorphic_multiplication(ct1, ct2);
    ct_mul_relin = Evaluator.relinearization_1(ct_mul_relin);
    poly_utils_1 mt_mul_relin = Evaluator.decryption(ct_mul_relin); // Use regular decryption for 2-element ciphertext
    int64_t nr_mul_relin = Evaluator.int_decode(mt_mul_relin);
    
    std::cout << "--- Performing ct_mul = Enc(m1) * Enc(m2) (with relinearization v1)" << std::endl;
    std::cout << "* ct_mul[0]: [";
    for (int64_t i = 0; i < 5; i++) {
        std::cout << ct_mul_relin[0].F[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dec = Dec(ct_mul)" << std::endl;
    std::cout << "* ct_dec: [";
    for (int64_t i = 0; i < 10; i++) {
        std::cout << mt_mul_relin.F[i];
        if (i < 9) std::cout << ", ";
    }
    std::cout << ", ...]" << std::endl;
    
    std::cout << "--- Performing ct_dcd = Decode(ct_dec)" << std::endl;
    std::cout << "* ct_dcd: " << nr_mul_relin << std::endl;
    
    if (nr_mul_relin == ne_mul) {
        std::cout << "* Homomorphic multiplication works." << std::endl;
    } else {
        std::cout << "* Homomorphic multiplication does not work." << std::endl;
    }
    std::cout << std::endl;
    
    return 0;
}