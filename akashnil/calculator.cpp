/**
Copyright 2018 Chia Network Inc
Modifications copyright (C) 2019 Akashnil Dutta

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
**/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <gmpxx.h>
#include <cassert>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <errno.h>

using namespace std;

#define LOG2(X) (63 - __builtin_clzll((X)))

inline uint64_t signed_shift(uint64_t op, int shift) {
    if (shift > 0) return op << shift;
    if (shift <= -64) return 0;
    return op >> (-shift);
}

uint64_t simple_hash(uint64_t x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x += 13962387127571770757ULL;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

struct form {
    mpz_t a;
    mpz_t b;
    mpz_t c;

    uint64_t hash() {
        return simple_hash(mpz_get_ui(a)) ^ simple_hash(mpz_get_ui(b));
    }
};

ostream& operator<<(ostream& os, const form& f) {
    return os << "a: " <<  f.a << "\n" << "b: " << f.b << "\n" << "c: " << f.c << "\n";
}

__int128 low_mask = 0xFFFFFFFFFFFFFFFF;

void mpz_set_128(mpz_t& rop, __int128 op) {
    bool neg = false;
    if(op < 0) {op = -op; neg = true;}
    mpz_set_ui(rop, op >> 64);
    mpz_mul_2exp(rop, rop, 64);
    mpz_add_ui(rop, rop, op & low_mask);
    if (neg) mpz_neg(rop, rop);
}

__int128 mpz_get_128(const mpz_t op) {
    auto size = mpz_size(op);
    // size must be <= 2
    __int128 hi = size > 1 ? mpz_getlimbn(op, 1) : 0;
    __int128 lo = mpz_getlimbn(op, 0); // size > 0 check?
    lo += hi << 64;
    return mpz_sgn(op) < 0 ? -lo : lo;
}

std::ostream&
operator<<( std::ostream& dest, __int128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}

form f_;

void normalize(form& f) {
    mpz_add(f_.b, f.b, f.c);
    mpz_mul_2exp(f_.c, f.c, 1);
    mpz_fdiv_qr (f_.b, f_.c, f_.b, f_.c);

    mpz_sub(f_.a, f.b, f.c);
    mpz_add(f_.a, f_.a, f_.c);

    mpz_sub(f.b, f.c, f_.c);

    mpz_div_2exp(f_.a, f_.a, 1);
    mpz_submul(f.a, f_.b, f_.a);

    mpz_swap(f.a, f.c);
}

void reduce(form& f) {
    if (mpz_cmp(f.a, f.c) < 0) {
        swap(f.a, f.c);
        mpz_neg (f.b, f.b);
    }
    while (mpz_cmp(f.a, f.c) > 0) normalize(f);
    if (mpz_cmp_ui(f.a, 1) == 0) mpz_set_ui(f.b, 1);
}

inline int64_t mpz_get_si_2exp (signed long int *exp, const mpz_t op) {
    uint64_t size = mpz_size(op);
    uint64_t last = mpz_getlimbn(op, size - 1);
    uint64_t ret;
    int lg2 = LOG2(last) + 1;
    *exp = lg2; ret = signed_shift(last, 63 - *exp);
    if (size > 1) {
        *exp += (size-1) * 64;
        uint64_t prev = mpz_getlimbn(op, size - 2);
        ret += signed_shift(prev, -1 - lg2);
    }
    if (mpz_sgn(op) < 0) return - ((int64_t)ret);
    return ret;
}

inline bool test_reduction(form& f) {
    int a_b = mpz_cmpabs(f.a, f.b);
    int c_b = mpz_cmpabs(f.c, f.b);

    if (a_b < 0 || c_b < 0) return false;

    int a_c = mpz_cmp(f.a, f.c);

    if (a_c > 0) {
        mpz_swap(f.a, f.c); mpz_neg(f.b, f.b);
    }

    if (a_c == 0 && mpz_sgn(f.b) < 0) {
        mpz_neg(f.b, f.b);
    }
    
    return true;
}

const int64_t THRESH = 1UL<<31;
const int64_t EXP_THRESH = 31;
mpz_t faa, fab, fac, fba, fbb, fbc, fca, fcb, fcc;

inline void fast_reduce(form& f) {

    int64_t u, v, w, x, u_, v_, w_, x_;
    int64_t delta, gamma, sgn;
    int64_t a, b, c, a_, b_, c_;
    int64_t aa, ab, ac, ba, bb, bc, ca, cb, cc;
    signed long int a_exp, b_exp, c_exp, max_exp, min_exp;

    while (!test_reduction(f)) {
        
        a = mpz_get_si_2exp(&a_exp, f.a);
        b = mpz_get_si_2exp(&b_exp, f.b);
        c = mpz_get_si_2exp(&c_exp, f.c);

        max_exp = a_exp;
        min_exp = a_exp;

        if (max_exp < b_exp) max_exp = b_exp;
        if (min_exp > b_exp) min_exp = b_exp;

        if (max_exp < c_exp) max_exp = c_exp;
        if (min_exp > c_exp) min_exp = c_exp;

        if (max_exp - min_exp > EXP_THRESH) {
            normalize(f); continue;
        }
        max_exp++; // for safety vs overflow

        a >>= (max_exp - a_exp);
        b >>= (max_exp - b_exp);
        c >>= (max_exp - c_exp);

        u_ = 1; v_ = 0; w_ = 0; x_ = 1;

        do {
            u = u_; v = v_; w = w_; x = x_;
            delta = b >= 0 ? (b+c) / (c<<1) : - (-b+c) / (c<<1);
            a_ = c;
            c_ = c * delta;
            b_ = -b + (c_ << 1);
            gamma = b - c_;
            c_ = a - delta * gamma; // may be overflow here which we'll handle later

            a = a_; b = b_; c = c_;

            u_ = v;
            v_ = -u + delta * v;
            w_ = x;
            x_ = -w + delta * x;

        } while ((abs(v_) | abs(x_)) <= THRESH && a > c && c > 0);

        if ((abs(v_) | abs(x_)) <= THRESH) {
            u = u_; v = v_; w = w_; x = x_;
        }

        aa = u * u; ab = u * w; ac = w * w;
        ba = u * v << 1; bb = u * x + v * w; bc = w * x << 1;
        ca = v * v; cb = v * x; cc = x * x;

        mpz_mul_si(faa, f.a, aa);
        mpz_mul_si(fab, f.b, ab);
        mpz_mul_si(fac, f.c, ac);

        mpz_mul_si(fba, f.a, ba);
        mpz_mul_si(fbb, f.b, bb);
        mpz_mul_si(fbc, f.c, bc);

        mpz_mul_si(fca, f.a, ca);
        mpz_mul_si(fcb, f.b, cb);
        mpz_mul_si(fcc, f.c, cc);

        mpz_add(f.a, faa, fab);
        mpz_add(f.a, f.a, fac);

        mpz_add(f.b, fba, fbb);
        mpz_add(f.b, f.b, fbc);

        mpz_add(f.c, fca, fcb);
        mpz_add(f.c, f.c, fcc);
    }
}

// Multiply by the fixed constant (2, 1, X) where X is inferred to match discriminant of f
void multiply2(form& f) {
    bool a_parity = mpz_odd_p(f.a) != 0;
    auto mod_4 = mpz_get_si(f.b) % 4;
    bool b_parity = mod_4 == 3 || mod_4 == -1; // is it 4k-1?
    bool c_parity = mpz_odd_p(f.c) != 0;

    if (b_parity ? a_parity : c_parity) {
        // assert((f.a + f.b + f.c) & 1 == 0);
        mpz_add(f.c, f.c, f.a);
        mpz_add(f.c, f.c, f.b);
        mpz_div_2exp(f.c, f.c, 1);
        mpz_mul_2exp(f.a, f.a, 1);
        mpz_add(f.b, f.b, f.a);
    } else if (b_parity) {
        // assert(f.a & 1 == 0);
        mpz_div_2exp(f.a, f.a, 1);
        mpz_mul_2exp(f.c, f.c, 1);
    } else {
        // assert(f.c & 1 == 0);
        mpz_mul_2exp(f.a, f.a, 1);
        mpz_div_2exp(f.c, f.c, 1);
    }
    reduce(f);
}

void square(form& f) {

    mpz_gcdext(f_.a, f_.b, f_.c, f.a, f.b);
    mpz_mul(f_.c, f.c, f_.c);
    mpz_fdiv_qr(f_.a, f_.c, f_.c, f.a);

    mpz_mul(f.c, f_.b, f.c);
    mpz_addmul(f.c, f_.a, f.b);
    mpz_addmul(f.c, f_.c, f_.c);

    // New b
    mpz_mul_2exp(f_.a, f.a, 1);
    mpz_submul(f.b, f_.a, f_.c);

    // New a
    mpz_mul(f.a, f.a, f.a);

    fast_reduce(f);
}

mpz_t x, y, m, u, v, n, beta, gmma, f1a, f2a, z;

// Assume that f1.a and f2.b are relatively prime
void multiply(form& prod, form& f1, form& f2, mpz_t& abs_disc) {
    
    // solve b1 + a1 z === b2 mod a2
    mpz_gcdext(m, x, y, f1.a, f2.a);
    // x is a1 inverse
    // z is x (b2 - b1) mod a2

    mpz_sub(z, f2.b, f1.b);
    mpz_div_2exp(z, z, 1);
    mpz_mul(z, z, x);
    mpz_mod(z, z, f2.a);
    
    mpz_mul(prod.a, f1.a, f2.a);
    mpz_mul_2exp(z, z, 1);
    mpz_mul(z, z, f1.a);
    mpz_add(prod.b, f1.b, z);
    
    mpz_mul(prod.c, prod.b, prod.b);
    mpz_add(prod.c, prod.c, abs_disc);
    mpz_mul_2exp(x, prod.a, 2);
    mpz_divexact(prod.c, prod.c, x);
    
    fast_reduce(prod);
}

void identity(form& f, mpz_t& abs_disc) {
    mpz_set_ui(f.a, 1);
    mpz_set_ui(f.b, 1);
    mpz_set(f.c, abs_disc);
    mpz_add_ui(f.c, f.c, 1);
    mpz_div_2exp(f.c, f.c, 2);
}

void generator(form& f, mpz_t& abs_disc) {
    mpz_set_ui(f.a, 2);
    mpz_set_ui(f.b, 1);
    mpz_set(f.c, abs_disc);
    mpz_add_ui(f.c, f.c, 1);
    mpz_div_2exp(f.c, f.c, 3);
}

void pow(form& f, __int128_t exp, mpz_t& abs_disc) {
    identity(f, abs_disc);
    for (int i = 127; i >= 0; i--) {
        if (exp >> i != 0) square(f);
        if ((exp >> i) & 1) multiply2(f);
    }
}

void generate_powers(mpz_t& abs_disc, __int128_t start, __int128_t step, __int128_t count) {
    form f, fs, p;
    mpz_inits(f.a, f.b, f.c, fs.a, fs.b, fs.c, p.a, p.b, p.c, NULL);
    pow(f, start, abs_disc);
    pow(fs, step, abs_disc);
    for (__int128_t i = 0; i < count; i++) {
        auto exp = start + step * i;
        cout << hex << f.hash() << " " << exp << endl;
        // Assuming that fs.a is prime here
        mpz_mod(p.a, f.a, fs.a);
        if (mpz_cmp_ui(p.a, 0) == 0) {
            pow(f, exp + step, abs_disc);
        } else {
            // fs.a must be relatively prime to p.a here for multiply to work correctly
            multiply(p, f, fs, abs_disc);
            mpz_set(f.a, p.a);
            mpz_set(f.b, p.b);
            mpz_set(f.c, p.c);
            // pow(p, exp + step, abs_disc);
            // assert(mpz_cmp(p.a, f.a) == 0);
        }
    }
}

#include <unordered_map>
unordered_map<uint64_t, uint64_t> exp_mapping;

void populate_hash_table(mpz_t& abs_disc, uint64_t count) {
    form f;
    mpz_inits(f.a, f.b, f.c, NULL);
    generator(f, abs_disc);
    for (uint64_t i = 0; i < count; i++) {
        auto key = f.hash();
        // Linear probing in case there's a hash collision
        while (exp_mapping.find(key) != exp_mapping.end()) key++;
        exp_mapping[key] = 1 + (i << 1);
        multiply2(f);
        multiply2(f);
    }
}

void measure_hash_table_capacity() {
    uint64_t cap = 1<<27;
    for (uint64_t j = 0; j < cap; j++) exp_mapping[j] = j;
    while (true) {
        uint64_t new_cap = cap * 4 / 3;
        for (uint64_t j = cap; j < new_cap; j++) exp_mapping[j] = j;
        cap = new_cap;
        cerr << cap << endl;
    }
}

#include <fstream>

mpz_t file_exp;

void read_files_for_collisions(mpz_t& abs_disc) {
    string filename;
    form f;
    mpz_inits(f.a, f.b, f.c, file_exp, NULL);
    while (true) {
        cin >> filename;
        std::ifstream infile(filename);
        uint64_t hsh;
        string exp_str;
        while (infile >> hex >> hsh >> exp_str) {
            while (exp_mapping.find(hsh) != exp_mapping.end()) {
                mpz_init_set_str(file_exp, exp_str.c_str(), 0);
                __int128_t exp = mpz_get_128(file_exp);
                uint64_t other = exp_mapping[hsh];
                pow(f, exp - other, abs_disc);
                if (mpz_cmp_ui(f.a, 1) == 0) {
                    cout << exp - other << endl;
                }
                pow(f, exp + other, abs_disc);
                if (mpz_cmp_ui(f.a, 1) == 0) {
                    cout << exp + other << endl;
                }
                hsh++;
            }
        }
        infile.close();
        remove(filename.c_str());
    }
}

#include <math.h>


__int128_t moduler_exponentiate(__int128_t a, int64_t n, __int128_t d) {
    if (n == 0) return 1;
    if ((n & 1) == 0) {
        __int128_t x = moduler_exponentiate(a, n/2, d);
        return (x * x) % d;
    } else {
        return (moduler_exponentiate(a, n-1, d) * a) % d;
    }
}
    
bool is_quadratic_residue(mpz_t& abs_disc, uint64_t prime) {
    if (prime == 2) return true;
    mpz_set_ui(x, prime);
    mpz_fdiv_r(x, abs_disc, x);
    return moduler_exponentiate(prime - mpz_get_128(x), (prime - 1) / 2, prime) == 1;
}

double find_approx(mpz_t& abs_disc, uint64_t start, uint64_t end) {
    bool* prime_sieve = (bool*) malloc(sizeof(bool) * (end - start));
    size_t base_size = sqrt(end) + 2;
    bool* base_sieve = (bool*) malloc(sizeof(bool) * base_size);
    for (uint64_t i = 0; i < base_size; i++) {
        base_sieve[i] = true;
    }
    for (uint64_t i = 0; i < end - start; i++) {
        prime_sieve[i] = true;
    }
    for (uint64_t i = 0; i < 2; i++) {
        base_sieve[i] = false;
        if (i >= start && i - start < end) prime_sieve[i - start] = false;
    }
    for (uint64_t i = 2; i < sqrt(base_size); i++) {
        if (!base_sieve[i]) continue;
        for (uint64_t j = i * i; j < base_size; j+= i) {
            base_sieve[j] = false;
        }
    }
    for (uint64_t i = 2; i < base_size; i++) {
        if (!base_sieve[i]) continue;
        auto u = start + i - 1; u -= u%i;
        if (u < i * i) u = i * i;
        for (uint64_t j = u; j < end; j += i) {
            prime_sieve[j - start] = false;
        }
    }
    double ret = 1.0;
    for (uint64_t i = 0; i < end - start; i++) {
        if (!prime_sieve[i]) continue;
        auto p = i + start;
        if (is_quadratic_residue(abs_disc, p)) ret += ret / p;
        else ret -= ret / p;
    }
    free (prime_sieve);
    free (base_sieve);
    return ret;
}

mpz_t abs_disc, start, step, num_steps;

/*
Usage:

1) ./a.out giant [discriminant] [start] [step] [num_steps]
Calculate g^(start + i * step) for i in range(num_steps). Output them to stdout.
g = (2,1,X). Assume that start is even, g^start = (a,b,c) has a prime.

2) ./a.out baby [discriminant] [num_steps]
Calculate h = g^(2i+1) for i in range(num_steps). Store (hash(h), 2i+1) in an unordered map.
Then keep reading for filenames. Each line in filename must be a pair (hash(g^exp), exp).
Try to find two hashes equal. In that case g^(exp1) = g^(exp2) or g^(exp1) = g^(exp2).
Then we must have g^(exp1 +- exp2) = identity. If such is found, print to stdout.

3) ./a.out approx [discriminant] [start] [end]
For all primes q in range (start, end) find out if -discriminant is a quadratic residue or not: 
F(q) = 2 or 0 respectively. Output one double number equal to prod_q (1 + (F(q)-1)/q)

4) ./a.out test_memory
Find out the largest number of elements that can be inserted into an unordered_map
uint64_t -> uint64_t. This is used to find the maximum value of [num_steps] in baby.
*/
int main(int argc, char* argv[]) {

    mpz_inits(f_.a, f_.b, f_.c, NULL);
    mpz_inits(x, y, m, u, v, n, beta, gmma, f1a, f2a, z, NULL);
    mpz_inits(faa, fab, fac, fba, fbb, fbc, fca, fcb, fcc, NULL);
    if (strcmp(argv[1], "giant") == 0) {
        mpz_inits(abs_disc, start, step, num_steps, NULL);
        mpz_init_set_str(abs_disc, argv[2], 0);
        mpz_init_set_str(start, argv[3], 0);
        mpz_init_set_str(step, argv[4], 0);
        mpz_init_set_str(num_steps, argv[5], 0);
        generate_powers(abs_disc, mpz_get_128(start), mpz_get_128(step), mpz_get_128(num_steps));
    }
    if (strcmp(argv[1], "baby") == 0) {
        mpz_init_set_str(abs_disc, argv[2], 0);
        mpz_init_set_str(num_steps, argv[3], 0);
        populate_hash_table(abs_disc, mpz_get_128(num_steps));
        cout << "populated" << endl;
        read_files_for_collisions(abs_disc);
    }

    if (strcmp(argv[1], "approx") == 0) {
        mpz_init_set_str(abs_disc, argv[2], 0);
        uint64_t start = stoull(argv[3]);
        uint64_t end = stoull(argv[4]);
        cout.precision(17);
        cout << find_approx(abs_disc, start, end);
    }
    if (strcmp(argv[1], "test_memory") == 0) {
        measure_hash_table_capacity();
    }

    return 0;
}
