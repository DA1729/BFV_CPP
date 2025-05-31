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
#include "bfv_lib.h"

using namespace std;

int64_t mod(int64_t a, int64_t m){
    int64_t r = a % m;
    return r < 0 ? r + m : r;
}

inline uint64_t round_to_int(double x) {
    return static_cast<uint64_t>(std::round(x));
}


bfv::bfv(size_t n, uint64_t q, uint64_t t, uint64_t T, uint64_t l, uint64_t p,
         double mu, double sigma, const std::vector<std::vector<uint64_t>>& qnp)
    : n(n), q(q), t(t), T(T), l(l), p(p), mu(mu), sigma(sigma), qnp(qnp),
      s_k(n, q, qnp)
{}

void bfv::secret_key_gen(){
    poly_tools s(n, q, qnp);
    s.randomize(2);
    s_k = s;
}

void bfv::public_key_gen(){
    poly_tools a(n, q, qnp);
    poly_tools e(n, q, qnp);

    a.randomize(q);
    e.randomize(0, false, 1, mu, sigma);

    poly_tools p_k_0 = -(a * s_k + e);
    poly_tools p_k_1 = a;

    p_k.clear();
    p_k.push_back(p_k_0);
    p_k.push_back(p_k_1);
}

void bfv::eval_key_gen_1(uint64_t T_in){
    T = T_in;
    l = static_cast<uint64_t>(std::floor(std::log(q)/std::log(T)));

    rl_k_1.clear();

    poly_tools sk_squared = s_k * s_k;

    for (uint64_t i = 0; i <= l; i++){
        poly_tools ai(n, q, qnp);
        poly_tools ei(n, q, qnp);

        ai.randomize(q);
        ei.randomize(0, false, 1, mu, sigma);

        poly_tools ts_2(n, q, qnp);

        ts_2.F.resize(n);

        for (size_t j = 0; j < n; j++){
            ts_2.F[j] = (mod_pow(T, i, q) * sk_squared.F[j]) % q;
        }

        poly_tools rl_ki0 = ts_2 - (ai * s_k + ei);
        poly_tools rl_ki1 = ai;

        rl_k_1.emplace_back(rl_ki0, rl_ki1);



    }
}


void bfv::eval_key_gen_2(uint64_t p_in){
    p = p_in;
    uint64_t pq = p * q;

    rl_k_2.clear();

    poly_tools a(n, pq, qnp);
    poly_tools e(n, pq, qnp);

    a.randomize(pq);
    e.randomize(0, false, 1, mu, sigma);


    std::vector<uint64_t> c0 = red_pol_mul_2(a.F, s_k.F);

    for (size_t i = 0; i < n; ++i) {
        c0[i] = (c0[i] + e.F[i]) % pq;
    }

    // Step 3: c1 = p * sk^2
    std::vector<uint64_t> c1 = red_pol_mul_2(s_k.F, s_k.F);
    for (size_t i = 0; i < n; ++i) {
        c1[i] = (p * c1[i]) % pq;
    }

    // Step 4: c2 = (c1 - c0) % pq
    std::vector<uint64_t> c2(n);
    for (size_t i = 0; i < n; ++i) {
        c2[i] = (c1[i] + pq - c0[i]) % pq;
    }

    // Step 5: store result in poly_tools
    poly_tools c(n, pq, qnp);
    c.F = c2;

    rl_k_2.emplace_back(c);
    rl_k_2.emplace_back(a);
}

std::vector<poly_tools> bfv::encryption(const poly_tools& m_in){

    int64_t delta = static_cast<int64_t>(std::floor(q / t));

    poly_tools u(n, q, qnp), e1(n, q, qnp), e2(n, q, qnp);
    u.randomize(2);
    e1.randomize(0, false, 1, mu, sigma);
    e2.randomize(0, false, 1, mu, sigma);

    poly_tools md(n, q, qnp);
    for (int i = 0; i < n; i++) {
        md.F[i] = (delta * m_in.F[i]) % q;
    }

    poly_tools c0 = p_k[0] * u + e1 + md;
    poly_tools c1 = p_k[1] * u + e2;

    return {c0, c1};

}

poly_tools bfv::decryption(const std::vector<poly_tools> &c_t_in){
    poly_tools m = c_t_in[1] * s_k + c_t_in[0];

    for (int i = 0; i < n; i++) {
        m.F[i] = round_to_int((t * m.F[i]) / static_cast<double>(q)) % t;
    }

    poly_tools mr(n, t, qnp);
    mr.F = m.F;
    mr.in_ntt = m.in_ntt;
    return mr;

    poly_tools sk2 = s_k * s_k;
    poly_tools m = c_t_in[0] + c_t_in[1] * s_k + c_t_in[2] * sk2;

    for (int i = 0; i < n; i++) {
        m.F[i] = round_to_int((t * m.F[i]) / static_cast<double>(q)) % t;
    }

    poly_tools mr(n, t, qnp);
    mr.F = m.F;
    mr.in_ntt = m.in_ntt;
    return mr;

}


poly_tools bfv::decryption_2(const std::vector<poly_tools> &c_t_in){
    poly_tools sk2 = s_k * s_k;
    poly_tools m = c_t_in[0] + c_t_in[1] * s_k + c_t_in[2] * sk2;

    for (int i = 0; i < n; i++) {
        m.F[i] = round_to_int((t * m.F[i]) / static_cast<double>(q)) % t;
    }

    poly_tools mr(n, t, qnp);
    mr.F = m.F;
    mr.in_ntt = m.in_ntt;
    return mr;
}

std::vector<poly_tools> bfv::relinearization_1(const std::vector<poly_tools> &ct){
    poly_tools c0 = ct[0];
    poly_tools c1 = ct[1];
    poly_tools c2 = ct[2];

    std::vector<poly_tools> c2i;
    poly_tools c2q(n, q, qnp);
    c2q.F = c2.F;

    for (int i = 0; i <= l; i++) {
        poly_tools c2r(n, q, qnp);
        for (int j = 0; j < n; j++) {
            int qt = c2q.F[j] / T;
            int rt = c2q.F[j] - qt * T;
            c2q.F[j] = qt;
            c2r.F[j] = rt;
        }
        c2i.push_back(c2r);
    }

    poly_tools c0r = c0;
    poly_tools c1r = c1;

    for (int i = 0; i <= l; i++) {
        c0r = c0r + (rl_k_1[i].first * c2i[i]);
        c1r = c1r + (rl_k_1[i].second * c2i[i]);
    }

    return {c0r, c1r};

}

std::vector<poly_tools> bfv::relinearization_1(const std::vector<poly_tools> &ct){
    poly_tools c0 = ct[0];
    poly_tools c1 = ct[1];
    poly_tools c2 = ct[2];

    std::vector<uint64_t> c2_0 = red_pol_mul_2(c2.F, rl_k_2[0].F);
    std::vector<uint64_t> c2_1 = red_pol_mul_2(c2.F, rl_k_2[1].F);

    for (auto &x : c2_0) x = round_to_int(x / static_cast<double>(p)) % q;
    for (auto &x : c2_1) x = round_to_int(x / static_cast<double>(p)) % q;

    poly_tools c0e(n, q, qnp); c0e.F = c2_0;
    poly_tools c1e(n, q, qnp); c1e.F = c2_1;

    poly_tools c0r = c0e + c0;
    poly_tools c1r = c1e + c1;

    return {c0r, c1r};
}

poly_tools bfv::int_encode(int64_t m){
    poly_tools mr(n, q, qnp);
    if (m > 0) {
        int64_t mt = m;
        for (int i = 0; i < n; i++) {
            mr.F[i] = mt % 2;
            mt /= 2;
        }
    } else if (m < 0) {
        int64_t mt = -m;
        for (int i = 0; i < n; i++) {
            mr.F[i] = (t - (mt % 2)) % t;
            mt /= 2;
        }
    }
    return mr;

}

int64_t bfv::int_decode(const poly_tools& m){
    int64_t result = 0;
    int64_t threshold = (t == 2) ? 2 : ((t + 1) >> 1);
    for (int i = 0; i < n; i++) {
        int64_t c = m.F[i];
        int64_t c_ = (c >= threshold) ? -(t - c) : c;
        result += (c_ * (1LL << i));
    }
    return result;
}

// homomorphic evaluations
std::vector<poly_tools> bfv::homomorphic_addition(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1){
    return {ct0[0] + ct1[0], ct0[1] + ct1[1]};
}

std::vector<poly_tools> bfv::homomorphic_subtraction(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1){
    return {ct0[0] - ct1[0], ct0[1] - ct1[1]};
}

std::vector<poly_tools> bfv::homomorphic_multiplication(const std::vector<poly_tools> &ct0, const std::vector<poly_tools> &ct1){
    std::vector<uint64_t> r0 = red_pol_mul_2(ct0[0].F, ct1[0].F);
    std::vector<uint64_t> r1 = red_pol_mul_2(ct0[0].F, ct1[1].F);
    std::vector<uint64_t> r2 = red_pol_mul_2(ct0[1].F, ct1[0].F);
    std::vector<uint64_t> r3 = red_pol_mul_2(ct0[1].F, ct1[1].F);

    std::vector<int64_t> c0(r0.size()), c1(r1.size()), c2(r3.size());
    for (size_t i = 0; i < r0.size(); ++i) c0[i] = round_to_int((t * r0[i]) / static_cast<double>(q)) % q;
    for (size_t i = 0; i < r1.size(); ++i) c1[i] = round_to_int((t * (r1[i] + r2[i])) / static_cast<double>(q)) % q;
    for (size_t i = 0; i < r3.size(); ++i) c2[i] = round_to_int((t * r3[i]) / static_cast<double>(q)) % q;

    poly_tools p0(n, q, qnp), p1(n, q, qnp), p2(n, q, qnp);
    p0.F = std::vector<uint64_t>(c0.begin(), c0.end());
    p1.F = std::vector<uint64_t>(c1.begin(), c1.end());;
    p2.F = std::vector<uint64_t>(c2.begin(), c2.end());

    return {p0, p1, p2};

}