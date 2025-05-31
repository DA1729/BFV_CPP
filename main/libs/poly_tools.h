// mainly dealing with n (degree + 1), 
//q (modulo), 
//ntt parameters (w, w_inv, psi, psi_inv)

#ifndef POLY_TOOLS_H
#define POLY_TOOLS_H
#include <vector>
#include <cstdint>
#include <cstddef>
#include <random>
#include <chrono>
#include <iostream>
#include "mod_tools.h"
#include "ntt_utils.h"
#include "prime_utils.h"



//
class poly_tools{
    public:
        poly_tools(size_t n, uint64_t q, const std::vector<std::vector<uint64_t>>& n_p);
//
        void randomize(int B, bool domain = false, int type = 0, double mu = 0, double sigma = 0);
        
        std::string to_string() const;

        
        //friend std::ostream& operator<<(std::ostream& os, const poly_tools& p);
        //std::ostream& operator<<(std::ostream& os, const poly_tools& p);
        poly_tools operator+(const poly_tools& b) const;
        poly_tools operator-(const poly_tools& b) const;
        poly_tools operator*(const poly_tools& b) const;
        poly_tools operator%(uint64_t base) const;
        poly_tools operator-() const;
        bool operator==(const poly_tools& b) const;

        poly_tools round() const;


        poly_tools to_ntt() const;
        poly_tools to_pol() const;
    
        size_t n;
        uint64_t q;
        std::vector<std::vector<uint64_t>> n_p;
        std::vector<uint64_t> F;
        bool in_ntt;

};

#endif