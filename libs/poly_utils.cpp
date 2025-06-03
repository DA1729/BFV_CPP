#include "poly_utils.h"
#include <random>
#include <cmath>
#include <sstream>

poly_utils_1 poly_init(int64_t n, int64_t q, const ntt_params& n_p){
    poly_utils_1 poly;
    poly.n = n;
    poly.q = q;
    poly.F = std::vector<int64_t>(n, 0);
    poly.n_p = n_p;
    poly.in_ntt = false;
    return poly;
}


void poly_randomize(poly_utils_1& poly, int64_t B, bool domain, int type, double mu, double sigma){
    std::random_device rd;
    std::mt19937 gen(rd());
    
    if (type == 0){
        std::uniform_int_distribution<int64_t> dist(-(B/2), B/2);
        for (size_t i = 0; i < poly.n; i++){
            poly.F[i] = ((dist(gen) % poly.q) + poly.q) % poly.q;
        }
    }
    else{
        std::normal_distribution<double> dist(mu, sigma);

        for (int64_t i = 0; i < poly.n; i++){
            int64_t val = static_cast<int64_t>(std::round(dist(gen)));
            poly.F[i] = ((val % poly.q) + poly.q) % poly.q;
        }
    }

    poly.in_ntt = domain;
}

std::string poly_to_string(const poly_utils_1& poly){
    std::ostringstream oss;
    int64_t tmp = std::min(poly.n, int64_t(8));

    oss << poly.F[0];
    for (int64_t i = 1; i < tmp; ++i) {
        oss << " + " << poly.F[i] << "*x^" << i;
    }

    if (poly.n > 8) {
        oss << " + ...";
    }

    return oss.str();
}

poly_utils_1 poly_add(const poly_utils_1& a, const poly_utils_1& b){
    if (a.in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial addition: inputs must be in the same domain");
    }

    else if(a.q != b.q){
        throw std::invalid_argument("polynomial addition: inputs must have the same modulus");
    }

    else{
        poly_utils_1 c = poly_init(a.n, a.q, a.n_p);
        for (size_t i = 0; i < a.n; i++){
            c.F[i] = (a.F[i] + b.F[i]) % a.q;
        }

        c.in_ntt = a.in_ntt;

        return c;
    }

    
}


poly_utils_1 poly_sub(const poly_utils_1& a, const poly_utils_1& b){
    if (a.in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial addition: inputs must be in the same domain");
    }

    else if(a.q != b.q){
        throw std::invalid_argument("polynomial addition: inputs must have the same modulus");
    }

    else{
        poly_utils_1 c = poly_init(a.n, a.q, a.n_p);

        for (size_t i = 0; i < a.n; i++){
            c.F[i] = (a.F[i] - b.F[i] + a.q) % a.q;
        }
        c.in_ntt = a.in_ntt;
        return c;
    }
}

poly_utils_1 poly_mul(const poly_utils_1& a, const poly_utils_1& b){
    if (a.in_ntt != b.in_ntt){
        throw std::invalid_argument("polynomial addition: inputs must be in the same domain");
    }

    else if(a.q != b.q){
        throw std::invalid_argument("polynomial addition: inputs must have the same modulus");
    }

    else{
        int64_t n = a.n;
        int64_t q = a.q;
        const std::vector<int64_t>& w = a.n_p.w;
        const std::vector<int64_t>& w_inv = a.n_p.w_inv;
        const std::vector<int64_t>& psi = a.n_p.psi;
        const std::vector<int64_t>& psi_inv = a.n_p.psi_inv;

        poly_utils_1 c = poly_init(n, q, a.n_p);

        if (a.in_ntt && b.in_ntt){
            for (size_t i = 0; i < n; i++){
                c.F[i] = (a.F[i]*b.F[i]) % q;
            }
            c.in_ntt = true;
            return c;
        }

        else{
            std::vector<int64_t> s_p(n), b_p(n), s_n(n), b_n(n), sb_n(n), sb_p(n), sb(n);

            for (size_t i = 0; i < n; i++){
                s_p[i] = (a.F[i] * psi[i]) % q;
                b_p[i] = (b.F[i] * psi[i]) % q;
            }

            s_n = ntt(s_p, w, q);
            b_n = ntt(b_p, w, q);

            for (size_t i = 0; i < n; i++){
                sb_n[i] = (s_n[i] * b_n[i]) % q;
            }

            sb_p = intt(sb_n, w_inv, q);

            for (size_t i = 0; i < n; i++){
                sb[i] = (sb_p[i] * psi_inv[i]) % q;
            }

            c.F = sb;
            c.in_ntt = false;

            return c;
        }
    }


}

poly_utils_1 poly_mod(const poly_utils_1& a, int64_t base){
    poly_utils_1 b = poly_init(a.n, a.q, a.n_p);

    for (size_t i = 0; i < a.n; i++){
        b.F[i] = a.F[i] % base;
    }

    b.in_ntt = a.in_ntt;

    return b;
}

poly_utils_1 poly_round(const poly_utils_1& a){
    poly_utils_1 b = poly_init(a.n, a.q, a.n_p);

    for (size_t i = 0; i < a.n; i++){
        b.F[i] = std::llround(a.F[i]);
    }

    b.in_ntt = a.in_ntt;

    return b;
}

bool poly_eq(const poly_utils_1& a, const poly_utils_1& b){
    if (a.n != b.n || a.q != b.q) return false;

    else{
        for (size_t i = 0; i < a.n; i++){
            if (a.F[i] != b.F[i]) return false;
        }

        return true;
    }
}

poly_utils_1 poly_neg(const poly_utils_1& a){
    poly_utils_1 b;
    b.n = a.n;
    b.q = a.q;
    b.F.resize(a.n);
    b.n_p = a.n_p;
    b.in_ntt = a.in_ntt;

    for (size_t i = 0; i < a.n; i++){
        b.F[i] = (a.F[i]) % a.q;
        if (b.F[i] < 0) b.F[i] += a.q;
    }

    return b;
}

poly_utils_1 poly_copy(poly_utils_1& dest, const poly_utils_1& a){
    dest.n = a.n;
    dest.q = a.q;
    dest.n_p = a.n_p;
    dest.in_ntt = a.in_ntt;
    dest.F = a.F;
    return dest;
}

poly_utils_1 to_ntt(const poly_utils_1& a){
    poly_utils_1 b; 
    b.n = a.n;
    b.q = a.q;
    b.n_p = a.n_p;
    b.F.resize(a.n);

    if (!a.in_ntt){
        b.F = ntt(a.F, a.n_p.w, a.q);
        b.in_ntt = true;
    }

    else {
        b.F = a.F;
        b.in_ntt = true;
    }

    return b;
}

poly_utils_1 to_pol(const poly_utils_1& a){
    poly_utils_1 b;
    b.n = a.n;
    b.q = a.q;
    b.n_p = a.n_p;
    b.F.resize(a.n);

    if (!a.in_ntt){
        b.F = a.F;
        b.in_ntt = false;
    }

    else{
        b.F = intt(a.F, a.n_p.w_inv, a.q);
        b.in_ntt = false;
    }

    return b;
}