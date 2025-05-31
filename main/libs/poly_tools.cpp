#include <cmath>
#include <stdexcept>
#include <ctime>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include "poly_tools.h"

poly_tools::poly_tools(size_t n, uint64_t q, const std::vector<std::vector<uint64_t>>& np)
    : n(n), q(q), n_p(n_p), F(n, 0), in_ntt(false) {}

void poly_tools::randomize(int B, bool domain = false, int type = 0, double mu = 0, double sigma = 0){
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    if (type == 0){
        std::uniform_int_distribution<int64_t> dist(-(B/2), B/2);
        for (size_t i = 0; i < n; i++){
            int64_t val = dist(rng);
            F[i] = (val % q + q)%q;
        }
    }
    else {
        std::normal_distribution<double> dist(mu, sigma);
        for (size_t i = 0; i < n; i++){
            int64_t val = static_cast<int64_t>(std::round(dist(rng)));
            F[i] = (val % q + q) % q;
        }
    }

    in_ntt = domain;
}

std::string poly_tools::to_string() const{
    size_t tmp = std::min(n, size_t(8));
    std::string res = std::to_string(F[0]);

    for (size_t i = 1; i < tmp; i++){
        res += " + " + std::to_string(F[i]) + "*x" + std::to_string(i);
    }

    if (n > 8) res += "+ ...";

    return res;

}

std::ostream& operator<<(std::ostream& os, const poly_tools& p){
    return os << p.to_string();
}


poly_tools poly_tools::operator+(const poly_tools& b) const{
    if (in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial addition: inputs must be in the same domain");

    }

    if (q != b.q){
        throw std::invalid_argument("polynomial addition: inputs must have the same modulus");
    }

    poly_tools result(n, q, n_p);
    for (size_t i = 0; i < n; i++){
        result.F[i] = (F[i] + b.F[i]) % q;

    }

    result.in_ntt = in_ntt;

    return result;
}


poly_tools poly_tools::operator-(const poly_tools& b) const{
    if (in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial subtraction: inupts must be in the same domain");
    }

    if (q != b.q){
        throw std::invalid_argument("polynomial subtraction: inputs must have the same modulus");
    }

    poly_tools result(n, q, n_p);

    for (size_t i = 0; i < n; i++){
        result.F[i] = (F[i] - b.F[i]) % q;
    }

    result.in_ntt = in_ntt;
    return result;
}

poly_tools poly_tools::operator*(const poly_tools& b) const{
    if (in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial multiplication: inputs must be in the same domain");
    }

    if (q != b.q){
        throw std::invalid_argument("polynomial multiplication: inputs must have the same modulus");
    }

    poly_tools result(n, q, n_p);
    
    // multiplication in ntt domain
    if (in_ntt == true && b.in_ntt == true){
        for (size_t i = 0; i < n; i++){
            result.F[i] = (F[i] * result.F[i]) % q;
            result.in_ntt = true;
        }
    }

    // else multiply in polynomial domain 
    else{

        // ntt parameters
        const std::vector<uint64_t>& w_table = n_p[0];
        const std::vector<uint64_t>& wv_table = n_p[1];
        const std::vector<uint64_t>& psi_table = n_p[2];
        const std::vector<uint64_t>& psi_inv_table = n_p[3];

        std::vector<uint64_t> s_p(n), b_p(n);
        for (size_t i = 0; i < n; i++){
            s_p[i] = (F[i] * psi_table[i]) % q;
            b_p[i] = (b.F[i] * psi_table[i]) % q;
        }

        std::vector<uint64_t> s_n = ntt(s_p, w_table, q);
        std::vector<uint64_t> b_n = ntt(b_p, w_table, q);

        std::vector<uint64_t> sb_n(n);

        for (size_t i = 0; i < n; i++){
            sb_n[i] = (s_n[i] * b_n[i]) % q;
        }

        std::vector<uint64_t> sb_p = intt(sb_n, wv_table, q);

        std::vector<uint64_t> sb(n);

        for (size_t i = 0; i < n; i++){
            sb[i] = (sb_p[i] * psi_inv_table[i]) % q;
        }

        result.F = sb;
        result.in_ntt = false;

    }

    return result;

}


poly_tools poly_tools::operator%(uint64_t base) const{
    poly_tools result(n, q, n_p);

    for (size_t i = 0; i < n; i++){
        result.F[i] = F[i] % base;
    }
    result.in_ntt = in_ntt;

    return result;
}

poly_tools poly_tools::round() const{
    poly_tools result(n, q, n_p);

    for (size_t i = 0; i < n; i++){
        result.F[i] = static_cast<uint64_t>(std::llround(static_cast<double>(F[i])));
    }

    result.in_ntt = in_ntt;

    return result;
}

bool poly_tools::operator==(const poly_tools& b) const{
    if (n != b.n){
        return false;
    }

    else if (q != b.q){
        return false;
    }

    else {
        for (size_t i = 0; i < n; i++){
            if (F[i] != b.F[i]){
                return false;
            }
        }

        return true;
    }
}

poly_tools poly_tools::operator-() const{
    poly_tools result(n, q, n_p);
    
    for (size_t i = 0; i < n; i++){
        result.F[i] = ((q-F[i]) % q);
    }

    result.in_ntt = in_ntt;

    return result;
}


poly_tools poly_tools::to_ntt() const{
    poly_tools result(n, q, n_p);

    if (in_ntt == false){
        result.F = ntt(F, n_p[0], q);
        result.in_ntt = true;
    }

    else{
        result.F = F;
        result.in_ntt = true;
    }

    return result;

}

poly_tools poly_tools::to_pol() const{
    poly_tools result(n, q, n_p);

    if (in_ntt == true){
        result.F = intt(F, n_p[1], q);
        result.in_ntt = false;
    }

    else{
        result.F = F;
        result.in_ntt = false;
    }

    return result;
}