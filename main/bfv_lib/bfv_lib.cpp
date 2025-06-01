#include "bfv_lib.h"
#include "./../libs/poly_utils.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <sstream>
#include <string>

bfv::bfv(int64_t n, int64_t q, int64_t t, double mu, double sigma, ntt_params qnp)
    : n(n), q(q), t(t), T(0), l(0), p(0), mu(mu), sigma(sigma), qnp(qnp)
    {
        s_k = poly_init(n, q, qnp);

        p_k.clear();
        rl_k1.clear();
        rl_k2.clear();
    }

std::string bfv::to_string() const{
    std::ostringstream oss;

    oss << "\n--- Parameters:\n";
    oss << "n     : " << n << "\n";
    oss << "q     : " << q << "\n";
    oss << "t     : " << t << "\n";
    oss << "T     : " << T << "\n";
    oss << "l     : " << l << "\n";
    oss << "p     : " << p << "\n";
    oss << "mu    : " << mu << "\n";
    oss << "sigma : " << sigma << "\n";
    return oss.str();

}

void bfv::print_params(){
    std::cout << to_string() << std::endl;
}

void bfv::secret_key_gen(){

    poly_utils_1 s_k = poly_init(n, q, qnp);
    poly_randomize(s_k, 2);

}

void bfv::public_key_gen(){
    poly_utils_1 a, e;
    a = poly_init(n, q, qnp);
    e = poly_init(n, q, qnp);

    poly_randomize(a, q);
    poly_randomize(e, 0, false, 1, mu, sigma);

    poly_utils_1 temp = poly_init(n, q, qnp);

    temp = poly_mul(a, s_k);
    temp = poly_add(temp, e);
    temp = poly_neg(temp);

    poly_utils_1 p_k0 = poly_init(n, q, qnp), p_k1 = poly_init(n, q, qnp);

    poly_copy(p_k0, temp);
    poly_copy(p_k1, a);

    p_k.clear();
    p_k.push_back(p_k0);
    p_k.push_back(p_k1);
    
}

void bfv::eval_key_gen_1(){
    T = 1 << 4;
    l = static_cast<int>(std::floor(std::log(q) / std::log(T)));

    poly_utils_1 s_sq = poly_init(n, q, qnp);
    s_sq = poly_mul(s_k, s_k);

    rl_k1.clear();

    for (size_t i = 0; i < n; i++){
        poly_utils_1 a = poly_init(n, q, qnp), e = poly_init(n, q, qnp), ts2 = poly_init(n, q, qnp), temp0 = poly_init(n, q, qnp), temp1 = poly_init(n, q, qnp);

        poly_randomize(a, q);
        poly_randomize(e, 0, false, 1, mu, sigma);

        for (size_t j = 0; j < n; j++){
            ts2.F[j] = (static_cast<int64_t>(std::pow(T, i)) * s_sq.F[j]) % q;
        }

        temp0 = poly_mul(a, s_k);
        temp0 = poly_add(temp0, e);
        temp0 = poly_sub(ts2, temp0);

        poly_copy(temp1, a);
        rl_k1.push_back({temp0, temp1});
    }
}

void bfv::eval_key_gen_2(){
    rl_k2.clear();

    poly_utils_1 a = poly_init(n, p*q, qnp), e = poly_init(n, p*q, qnp), c = poly_init(n, p*q, qnp), temp1 = poly_init(n, p*q, qnp);

    poly_randomize(a, p*q);
    poly_randomize(e, 0, false, 1, mu, sigma);

    poly_utils_1 sk_sq = poly_init(n, p*q, qnp);
    sk_sq = poly_mul(s_k, s_k);

    temp1 = poly_mul(a, s_k);
    temp1 = poly_add(temp1, e);

    for (size_t i = 0; i < n; i++){
        temp1.F[i] = ((p * q) - (temp1.F[i] % (p * q))) % (p * q);
    }

    for (size_t i = 0; i < n; i++){
        c.F[i] = (temp1.F[i] + (p * sk_sq.F[i])) % (p * q);
    }

    rl_k2.push_back(c);
    rl_k2.push_back(a);

}


std::vector<poly_utils_1> bfv::encryption(const poly_utils_1& m){
    std::vector<poly_utils_1> cipher_text;

    int64_t delta = q/t;

    poly_utils_1 u = poly_init(n, q, qnp), e1 = poly_init(n, q, qnp), e2 = poly_init(n, q, qnp), md = poly_init(n, q, qnp), c0 = poly_init(n, q, qnp), c1 = poly_init(n, q, qnp);

    poly_randomize(u, 2);
    poly_randomize(e1, 0, false, 1, mu, sigma);
    poly_randomize(e2, 0, false, 1, mu, sigma);

    for (size_t i = 0; i < n; i++){
        md.F[i] = (delta * m.F[i]) % q;
    }

    c0 = poly_mul(p_k[0], u);

    for (size_t i = 0; i < n; i++){
        c0.F[i] = (c0.F[i] + e1.F[i] + md.F[i]) % q;
    }

    c1 = poly_mul(p_k[1], u);

    for (size_t i = 0; i < n; i++){
        c1.F[i] = (c1.F[i] + e2.F[i]) % q;
    }

    cipher_text.push_back(c0);
    cipher_text.push_back(c1);

    return cipher_text;
}

poly_utils_1 bfv::decryption(const std::vector<poly_utils_1>& c_t){
    poly_utils_1 m = poly_init(n, q, qnp), temp = poly_init(n, q, qnp);

    temp = poly_mul(c_t[1], s_k);
    m = poly_add(temp, c_t[0]);


    for (size_t i = 0; i < n; i++){
        double scaled = (double)(m.F[i]) * ((double)t / (double)q);
        m.F[i] = (int)std::round(scaled)%t;
        
        if (m.F[i] < 0) m.F[i] += t;
    }

    m.in_ntt = false;
    return m;


}

poly_utils_1 bfv::decryption_2(const std::vector<poly_utils_1>& c_t){
    poly_utils_1 sk2 = poly_mul(s_k, s_k);
    poly_utils_1 m = poly_add(c_t[0], poly_add(poly_mul(c_t[1], s_k), poly_mul(c_t[2], sk2)));

    for (size_t i = 0; i < n; i++){
        double scaled = ((double)t * m.F[i]) / q;
        m.F[i] = ((int64_t)std::round(scaled)) % t;
        if (m.F[i] < 0) m.F[i] += t;
    }

    poly_utils_1 mr = poly_init(n, t, qnp);
    mr.F = m.F;
    mr.in_ntt = m.in_ntt;

    return mr;  
}

std::vector<poly_utils_1> bfv::homomorphic_addition(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1){
    std::vector<poly_utils_1> result;
    poly_utils_1 c0 = poly_init(n, q, qnp), c1 = poly_init(n, q, qnp);

    for (size_t i = 0; i < n; i++){
        c0.F[i] = (c_t0[0].F[i] + c_t1[0].F[i]) % q;
        c1.F[i] = (c_t0[1].F[i] + c_t1[1].F[i]) % q;
    }

    result.push_back(c0);
    result.push_back(c1);

    return result;
}

std::vector<poly_utils_1> bfv::homomorphic_subtraction(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1){
    std::vector<poly_utils_1> result;
    poly_utils_1 c0 = poly_init(n, q, qnp), c1 = poly_init(n, q, qnp);

    for (size_t i = 0; i < n; i++){
        c0.F[i] = (c_t0[0].F[i] - c_t1[0].F[i] + q) % q;
        c1.F[i] = (c_t0[1].F[i] - c_t1[1].F[i] + q) % q;

    }

    result.push_back(c0);
    result.push_back(c1);

    return result;
}

std::vector<poly_utils_1> bfv::homomorphic_multiplication(const std::vector<poly_utils_1>& c_t0, const std::vector<poly_utils_1>& c_t1){
    poly_utils_1 c0 = poly_init(n, q, qnp), c1 = poly_init(n, q, qnp), c2 = poly_init(n, q, qnp);

    std::vector<int64_t> r0 = red_pol_mul_2(c_t0[0].F, c_t1[0].F);
    std::vector<int64_t> r1 = red_pol_mul_2(c_t0[0].F, c_t1[1].F);
    std::vector<int64_t> r2 = red_pol_mul_2(c_t0[1].F, c_t1[0].F);
    std::vector<int64_t> r3 = red_pol_mul_2(c_t0[1].F, c_t1[1].F);

    for (size_t i = 0; i < n; i++){
        double val_c0 = ((double)r0[i]) * ((double)t / (double)q);
        double val_c1 = ((double)(r1[i] + r2[i])) * ((double)t / (double)q);
        double val_c2 = ((double)r3[i]) * ((double)t / (double)q);

        c0.F[i] = (int64_t)std::round(val_c0) % q;
        c1.F[i] = (int64_t)std::round(val_c1) % q;
        c2.F[i] = (int64_t)std::round(val_c2) % q;

        if (c0.F[i] < 0) c0.F[i] += q;
        if (c1.F[i] < 0) c1.F[i] += q;
        if (c2.F[i] < 0) c2.F[i] += q;

    }

    return {c0, c1, c2};


}

// encoding and decoding
poly_utils_1 bfv::int_encode(int64_t m){
    poly_utils_1 encoded = poly_init(n, t, qnp);

    if (m > 0) {
        int mt = m;
        for (size_t i = 0; i < n; i++){
            encoded.F[i] = mt % 2;
            mt /= 2;
        }
    } else if (m < 0){
        int mt = -m;
        for (size_t i = 0; i < n; i++){
            encoded.F[i] = (t - (mt % 2)) % t;
            mt /= 2;
        }
    } else{

    }

    return encoded;
}

int64_t bfv::int_decode(const poly_utils_1& m){
    int result = 0;
    int64_t threshold = (t == 2) ? 2 : ((t + 1) >> 1);

    for (size_t i = 0; i < n; i++){
        int64_t c = m.F[i];
        int64_t c_ = (c >= threshold) ? -(t - c) : c;
        result += c_ * (1LL << i);
    }

    return result;
}


std::vector<poly_utils_1> bfv::relinearization_1(const std::vector<poly_utils_1>& c_t){
    poly_utils_1 c0 = c_t[0];
    poly_utils_1 c1 = c_t[1];
    poly_utils_1 c2 = c_t[2];

    std::vector<poly_utils_1> c2i;
    poly_utils_1 c2q = poly_init(n, q, qnp);
    c2q.F = c2.F;

    for (size_t i = 0; i <= l; i++){
        poly_utils_1 c2r = poly_init(n, q, qnp);
        for (size_t j = 0; j < n; j++){
            int64_t qt = c2q.F[j] / T;
            int64_t rt = c2q.F[j] - qt * T;

            c2q.F[j] = qt;
            c2r.F[j] = rt;
        }

        c2i.push_back(c2r);
    }

    poly_utils_1 c0r = poly_init(n, q, qnp), c1r = poly_init(n, q, qnp);
    c0r.F = c0.F;
    c1r.F = c1.F;

    for (size_t i = 0; i < l; i++){
        c0r = poly_add(c0r, poly_mul(rl_k1[i][0], c2i[i]));
        c1r = poly_add(c1r, poly_mul(rl_k1[i][1], c2i[i]));
        
    }

    return {c0r, c1r};
}

std::vector<poly_utils_1> bfv::relinearization_2(const std::vector<poly_utils_1>& c_t){
    poly_utils_1 c0 = c_t[0];
    poly_utils_1 c1 = c_t[1];
    poly_utils_1 c2 = c_t[2];

    std::vector<int64_t> c2_0 = red_pol_mul_2(c2.F, rl_k2[0].F);
    std::vector<int64_t> c2_1 = red_pol_mul_2(c2.F, rl_k2[1].F);

    for (size_t i = 0; i < n; i++){
        c2_0[i] = ((int64_t)std::round((double)c2_0[i] / p)) % q;
        c2_1[i] = ((int64_t)std::round((double)c2_1[i] / p)) % q;
    }

    poly_utils_1 c0e = poly_init(n, q, qnp), c1e = poly_init(n, q, qnp);
    c0e.F = c2_0;
    c1e.F = c2_1;                  
    
    poly_utils_1 c0r = poly_add(c0, c0e);
    poly_utils_1 c1r = poly_add(c1, c1e);

    return {c0r, c1r};
}