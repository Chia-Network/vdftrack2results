/**
Copyright 2018 Chia Network Inc

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
#include <gmpxx.h>
#include <cassert>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <bitset>
#include <array>
#include <algorithm>

#include <pthread.h>

using namespace std;

long TEN9 = 1000000000;

struct form {
	// y = ax^2 + bxy + y^2
	mpz_t a;
	mpz_t b;
	mpz_t c;

	mpz_t d; // discriminant
};

form x;
form prin;
form giant_step2;
mpz_t Q_final;
mpz_t baby2;

mpz_t L, D;
mpz_t discriminant;

const long babyn = 2000000;
long babytot = 0; //total amount of baby exponents covered thus far.
long babystart = 0; //starting position to dump additional baby entries.

//vector<pair<array<char, 25>, long>> babytable(babyn);
vector<long> babyltable(babyn);
vector<long> primetable(3000000);

struct gianttask {
	mpz_t iter;
	mpz_t upper_bound;
	form xp;
};

struct babytask {
	long lo;
	long hi;
};

struct primetask {
	long lo;
	long hi;
	mpf_t res;
};

gianttask gtasks[6];
babytask btasks[6];
primetask ptasks[6];
pthread_t slaves[6];

ostream& operator<<(ostream& os, const form& f) {
	return os << "a: " <<  f.a << endl << "b: " << f.b << endl << "c: " << f.c << endl;
}

inline void normalize(form& f) {
	mpz_t negative_a, denom, old_b, r, ra;
	mpz_inits(negative_a, denom, old_b, r, ra, NULL);

	mpz_neg(negative_a, f.a);
	if (mpz_cmp(f.b, negative_a) > 0 && mpz_cmp(f.b, f.a) <= 0) {
		// Already normalized
		mpz_clears(negative_a, denom, old_b, r, ra, NULL);
		return;
	}
	// r = (a - b) / 2a
	// a = a
	// b = b + 2ra
	// c = ar^2 + br + c
	mpz_sub(r, f.a, f.b);

	mpz_mul_ui(denom, f.a, 2);

	// r = (a-b) / 2a
	mpz_fdiv_q(r, r, denom);

	mpz_set(old_b, f.b);

	mpz_mul(ra, r, f.a);
	mpz_add(f.b, f.b, ra);
	mpz_add(f.b, f.b, ra);

	// c += ar^2
	mpz_mul(ra, ra, r);
	mpz_add(f.c, f.c, ra);

	// c += rb
	mpz_set(ra, r);
	mpz_mul(ra, ra, old_b);
	mpz_add(f.c, f.c, ra);

	mpz_clears(negative_a, denom, old_b, r, ra, NULL);
}

inline void reduce(form& f) {
	mpz_t s, x, old_a, old_b;
	mpz_inits(s, x, old_a, old_b, NULL);

	normalize(f);
	while ((mpz_cmp(f.a, f.c) > 0) ||
		   (mpz_cmp(f.a, f.c) == 0 && mpz_cmp_si(f.b, 0) < 0)) {
		mpz_add(s, f.c, f.b);

		// x = 2c
		mpz_mul_ui(x, f.c, 2);
		mpz_fdiv_q(s, s, x);

		mpz_set(old_a, f.a);
		mpz_set(old_b, f.b);

		// b = -b
		mpz_set(f.a, f.c);
		mpz_neg(f.b, f.b);

		// x = 2sc
		mpz_mul(x, s, f.c);
		mpz_mul_ui(x, x, 2);

		// b += 2sc
		mpz_add(f.b, f.b, x);

		// c = cs^2
		mpz_mul(f.c, f.c, s);
		mpz_mul(f.c, f.c, s);

		// x = bs
		mpz_mul(x, old_b, s);

		// c -= bs
		mpz_sub(f.c, f.c, x);

		// c += a
		mpz_add(f.c, f.c, old_a);
	}
	normalize(f);

	mpz_clears(s, x, old_a, old_b, NULL);
}

inline form generator_for_discriminant(mpz_t* d) {
	mpz_t denom;
	mpz_init(denom);

	form x;
	mpz_init_set_ui(x.a, 2);
	mpz_init_set_ui(x.b, 1);
	mpz_init(x.c);
	mpz_init_set(x.d, *d);

	// c = b*b - d
	mpz_mul(x.c, x.b, x.b);
	mpz_sub(x.c, x.c, x.d);

	// denom = 4a
	mpz_mul_ui(denom, x.a, 4);

	mpz_fdiv_q(x.c, x.c, denom);
	reduce(x);

	mpz_clear(denom);
	return x;
}

// only works because d is odd.
inline form generate_principal_form(mpz_t* d) {
	form x;
	mpz_init_set_ui(x.a, 1);
	mpz_init_set_ui(x.b, 1);
	mpz_init_set_ui(x.c, 1);
	mpz_init_set(x.d, *d);
	mpz_sub(x.c, x.c, x.d);
	mpz_divexact_ui(x.c, x.c, 4);
	return x;
}

// Returns mu and v, solving for x:  ax = b mod m
// such that x = u + vn (n are all integers). Assumes that mu and v are initialized.
// Returns 0 on success, -1 on failure
inline int solve_linear_congruence(mpz_t& mu, mpz_t& v, mpz_t& a, mpz_t& b, mpz_t& m) {
	// g = gcd(a, m), and da + em = g
	mpz_t g, d, e, q, r;
	mpz_inits(g, d, e, q, r, NULL);

	mpz_gcdext(g, d, e, a, m);

	// q = b/g, r = b % g
	mpz_fdiv_qr(q, r, b, g);

	if (mpz_cmp_ui(r, 0) != 0) {
		// No solution, return error. Optimize out for speed..
		cout << "No solution to congruence" << endl;
		mpz_clears(g, d, e, q, r, NULL);
		return -1;
	}

	mpz_mul(mu, q, d);
	mpz_mod(mu, mu, m);

	mpz_fdiv_q(v, m, g);

	mpz_clears(g, d, e, q, r, NULL);
	return 0;
}

// Faster version without check, and without returning v
inline int solve_linear_congruence(mpz_t& mu, mpz_t& a, mpz_t& b, mpz_t& m) {
	mpz_t g, d, e, q;
	mpz_inits(g, d, e, q, NULL);

	mpz_gcdext(g, d, e, a, m);
	mpz_fdiv_q(q, b, g);
	mpz_mul(mu, q, d);
	mpz_mod(mu, mu, m);

	mpz_clears(g, d, e, q, NULL);
	return 0;
}

// Takes the gcd of three numbers
inline void three_gcd(mpz_t& ret, mpz_t& a, mpz_t& b, mpz_t& c) {
	mpz_gcd(ret, a, b);
	mpz_gcd(ret, ret, c);
}

// f1 gets multiplied by f2.
inline void multiplyinplace(form &f1, form &f2) {
	mpz_t g, h, w, j, r, s, t, u, a, b, m, mu, v, lambda, sigma, k, l;
	mpz_inits(g, h, w, j, r, s, t, u, a, b, m, mu, v, lambda, sigma, k, l, NULL);

	assert(mpz_cmp(f1.d, f2.d) == 0);

	// g = (b1 + b2) / 2
	mpz_add(g, f1.b, f2.b);
	mpz_fdiv_q_ui(g, g, 2);


	// h = (b2 - b1) / 2
	mpz_sub(h, f2.b, f1.b);
	mpz_fdiv_q_ui(h, h, 2);

	// w = gcd(a1, a2, g)
	three_gcd(w, f1.a, f2.a, g);

	// j = w
	mpz_set(j, w);

	// r = 0
	mpz_set_ui(r, 0);

	// s = a1/w
	mpz_fdiv_q(s, f1.a, w);

	// t = a2/w
	mpz_fdiv_q(t, f2.a, w);

	// u = g/w
	mpz_fdiv_q(u, g, w);

	// solve (tu)k = (hu + sc1) mod st, of the form k = mu + vn

	// a = tu
	mpz_mul(a, t, u);

	// b = hu + sc1
	mpz_mul(b, h, u);
	mpz_mul(m, s, f1.c);
	mpz_add(b, b, m);

	// m = st
	mpz_mul(m, s, t);

	int ret = solve_linear_congruence(mu, v, a, b, m);

	assert(ret == 0);



	// solve (tv)n = (h - t * mu) mod s, of the form n = lamda + sigma n'

	// a = tv
	mpz_mul(a, t, v);

	// b = h - t * mu
	mpz_mul(m, t, mu); // use m as a temp variable
	mpz_sub(b, h, m);

	// m = s
	mpz_set(m, s);

	ret = solve_linear_congruence(lambda, sigma, a, b, m);
	assert(ret == 0);

	// k = mu + v*lamda
	mpz_mul(a, v, lambda); // use a as a temp variable

	mpz_add(k, mu, a);

	// l = (k*t - h) / s
	mpz_mul(l, k, t);
	mpz_sub(l, l, h);
	mpz_fdiv_q(l, l, s);

	// m = (tuk - hu - cs) / st
	mpz_mul(m, t, u);
	mpz_mul(m, m, k);
	mpz_mul(a, h, u); // use a as a temp variable
	mpz_sub(m, m, a);
	mpz_mul(a, f1.c, s); // use a as a temp variable
	mpz_sub(m, m, a);
	mpz_mul(a, s, t); // use a as a temp variable
	mpz_fdiv_q(m, m, a);

	// A = st - ru
	mpz_mul(f1.a, s, t);
	mpz_mul(a, r, u); // use a as a temp variable
	mpz_sub(f1.a, f1.a, a);

	// B = ju + mr - (kt + ls)
	mpz_mul(f1.b, j, u);
	mpz_mul(a, m, r); //use a as a temp variable
	mpz_add(f1.b, f1.b, a);
	mpz_mul(a, k, t); // use a as a temp variable
	mpz_sub(f1.b, f1.b, a);
	mpz_mul(a, l, s); // use a as a temp variable
	mpz_sub(f1.b, f1.b, a);

	// C = kl - jm
	mpz_mul(f1.c, k, l);
	mpz_mul(a, j, m);
	mpz_sub(f1.c, f1.c, a);

	reduce(f1);

	mpz_clears(g, h, w, j, r, s, t, u, a, b, m, mu, v, lambda, sigma, k, l, NULL);
}

inline void squareinplace(form &f1) {
	mpz_t mu, m, old_a, a;
	mpz_inits(mu, m, old_a, a, NULL);

	int ret = solve_linear_congruence(mu, f1.b, f1.c, f1.a);
	assert(ret == 0);

	mpz_mul(m, f1.b, mu);
	mpz_sub(m, m, f1.c);
	mpz_fdiv_q(m, m, f1.a);

	// New a
	mpz_set(old_a, f1.a);
	mpz_mul(f1.a, f1.a, f1.a);

	// New b
	mpz_mul(a, mu, old_a);
	mpz_mul_ui(a, a, 2);
	mpz_sub(f1.b, f1.b, a);

	// New c
	mpz_mul(f1.c, mu, mu);
	mpz_sub(f1.c, f1.c, m);
	//mpz_set(f1.d, f1.d);
	reduce(f1);

	mpz_clears(mu, m, old_a, a, NULL);
}

inline void exponentiateinplace(form& f, mpz_t& E) {
	mpz_t e;
	mpz_init_set(e, E);

	mpz_sub_ui(e, e, 1);

	form f2;
	mpz_inits(f2.a, f2.b, f2.c, f2.d, NULL);
	mpz_set(f2.a, f.a);
	mpz_set(f2.b, f.b);
	mpz_set(f2.c, f.c);
	mpz_set(f2.d, f.d);

	while(mpz_cmp_si(e, 0)) {
		if(!mpz_divisible_2exp_p(e, 1)) {
			multiplyinplace(f, f2);
		}
		squareinplace(f2);
		mpz_tdiv_q_2exp(e, e, 1);
	}

	mpz_clears(f2.a, f2.b, f2.c, f2.d, NULL);
	mpz_clear(e);
}

int max_baby_arr_fill = 0;
long tot_baby_arr_fill = 0;

long baby_encodelong(form& f) {
	return mpz_get_ui(f.a);
}

inline form form_copy(form& f) {
	form ret;
	mpz_init_set(ret.a, f.a);
	mpz_init_set(ret.b, f.b);
	mpz_init_set(ret.c, f.c);
	mpz_init_set(ret.d, f.d);
	return ret;
}

// NOT about being in same equivalance class.
inline bool form_equal(form& x, form& y){
	return !mpz_cmp(x.a, y.a) && !mpz_cmp(x.b, y.b) && !mpz_cmp(x.c, y.c) && !mpz_cmp(x.d, y.d);
}

char temp[1234567];

timespec timer;
long curtime = 0;
long newtime = 0;


void timestamp(){
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timer);
	newtime = TEN9 * timer.tv_sec + timer.tv_nsec;

	long pnewtime = newtime/6; //6 cores. alas its only like x4 but whatever.
	long pcurtime = curtime/6; //6 cores

	long seconds = (pnewtime - pcurtime)/1e9;
	long minutes = seconds / 60 % 60;
	long hours = seconds / 3600 % 3600;
	seconds = seconds % 60;

	cout << hours << " hours, " << minutes << " minutes, " << seconds << " seconds" << endl;
	curtime = newtime;
}


void* work(void* task_ptr) {

	gianttask* task_pointer = (gianttask*) task_ptr;

	mpz_t iter, upper_bound;
	mpz_inits(iter, upper_bound, NULL);
	mpz_set(iter, task_pointer->iter);
	mpz_set(upper_bound, task_pointer->upper_bound);
	form xp = form_copy(task_pointer->xp);

	for (; mpz_cmp(iter, upper_bound) <= 0; mpz_add(iter, iter, baby2)) {
		//array<char, 25> xps = baby_encodearr(xp);
		long xps = baby_encodelong(xp);

		auto match = lower_bound(babyltable.begin(), babyltable.end(), xps);

		if (match != babyltable.end() && xps == *match) {

			cout << "HIT (F/P) NEAR ORDER: " << iter << endl;
		}

		multiplyinplace(xp, giant_step2);
		reduce(xp);

	}

	return NULL;
}

void* babywork(void* task_ptr) {

	babytask* task_pointer = (babytask*) task_ptr;
	long lo = task_pointer->lo;
	long hi = task_pointer->hi;
	long loptr = lo + babystart;
	long hiptr = hi + babystart;
	long loexpbegin = lo + babytot; 

	form babyx;
	if (loexpbegin == 0) {
		babyx = form_copy(prin);
	} else {
		babyx = form_copy(x);
		mpz_t E;
		mpz_init_set_ui(E, loexpbegin);
		exponentiateinplace(babyx, E);
		mpz_clear(E);
	}

	for (long i=loptr; i<hiptr; i++) {
		//array<char, 25> babys = baby_encodearr(babyx);
		long babys = baby_encodelong(babyx);
		babyltable[i] = babys;
		multiplyinplace(babyx, x);
	}

	mpz_clears(babyx.a, babyx.b, babyx.c, babyx.d, NULL);

	return NULL;
}

void* primework(void* task_ptr) {

	primetask* task_pointer = (primetask*) task_ptr;
	long lo = task_pointer->lo;
	long hi = task_pointer->hi;

	mpz_t dis, pmod;
	mpf_t fone, fp, fleg, accum;
	mpz_inits(dis, pmod, NULL);
	mpf_inits(fone, fp, fleg, accum, NULL);
	mpf_set_ui(fone, 1);

	mpf_set_ui(task_pointer->res, 1);

	for (long i=lo; i<hi; i++) {

		long p = primetable[i];
		
		mpz_mod_ui(dis, discriminant, p);

		long leggg = 1;

		mpz_set_ui(pmod, p);

		if (p != 2) {
			leggg = mpz_legendre(dis, pmod);
		}

		mpf_set_ui(fp, p);
		mpf_set_si(fleg, leggg);
		mpf_div(accum, fleg, fp);
		mpf_neg(accum, accum);
		mpf_add_ui(accum, accum, 1);
		mpf_div(accum, fone, accum);
		mpf_mul(task_pointer->res, task_pointer->res, accum);
		//Q *= 1.0/(1.0-( (leg_symbol*1.0) / (p*1.0) ));
	}


	mpz_clears(dis, pmod, NULL);
	mpf_clears(fone, fp, fleg, accum, NULL);

	return NULL;
}

//USAGE == ./a.out {discriminant} {hexseed}
int main(int argc, char* argv[]) {

	//mpf_set_default_prec(64); //1 limbs
	mpf_set_default_prec(128); //2 limbs


	mpz_inits(L, D, NULL);

	mpf_t P, dabs, dsqrt, dsqrtoverpi, PIE;
	mpf_inits(P, dabs, dsqrt, dsqrtoverpi, PIE, NULL);

	// NOTE: UPDATE THIS WHEN YOU CHANGE PRIMEFILE.
	//string P_str = "1048576e0"; //2^20
	string P_str = "33554432e0"; //2^25
	//string P_str = "1073741824e0"; //2^30
	//string P_str = "2147483648e0"; //2^31
	//string P_str = "4294967296e0"; //2^32 
	//string P_str = "8589934592e0"; //2^33 
	//string P_str = "17179869184e0"; //2^34 
	//string P_str = "34359738368e0"; //2^35 
	//string P_str = "68719476736e0"; //2^36
	//string P_str = "137438953472e0"; //2^37
	//string P_str = "274877906944e0"; //2^38

	strcpy(temp, P_str.c_str());
	mpf_set_str(P, temp, 10);

	string pies = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862e0";
	strcpy(temp, pies.c_str());
	mpf_set_str(PIE, temp, 10);

	mpz_t dis;
	mpz_init(dis);

	assert(mpz_init_set_str(discriminant, argv[1], 0) == 0);

	mpz_set(D, discriminant);
	mpz_abs(L, D);
	mpz_root(L, L, 4);


	string hexseed(argv[2]);
	cout << "hexseed: " << hexseed << endl;

	//using generator_for_discriminant to create a non principal form is fine.
	x = generator_for_discriminant(&discriminant);
	cout << "starting form: " << x.a << " " << x.b << " " << x.c << endl;
	prin = generate_principal_form(&discriminant);
	
	stringstream prinb;
	prinb << prin;
	string prins = prinb.str();


	long REMAINING_BABY_DUP_REMOVES = 3;
	
	while (babystart != babyn) {
		long babytask_distance = (babyn-babystart)/6;

		cout << "babytot: " << babytot << " babystart: " << babystart << endl; 


		for (int i=0; i<6; i++) {
			btasks[i].lo = babytask_distance * i;
			if (i != 6-1) {
				btasks[i].hi = babytask_distance * (i+1);
			} else {
				btasks[i].hi = (babyn-babystart);
			}
		}

		for (int i=0; i<6; i++) {
			pthread_create(&slaves[i], NULL, babywork, &btasks[i]);
		}
		for (int i=0; i<6; i++) {
			pthread_join(slaves[i], NULL);
		}

		sort(babyltable.begin(), babyltable.end());

		babytot += (babyn - babystart);

		if (REMAINING_BABY_DUP_REMOVES == 0) break;

		babystart = (long)(unique(babyltable.begin(), babyltable.end()) - babyltable.begin());

		REMAINING_BABY_DUP_REMOVES--;
	}

	cout << "babytablesize (removed duplicate entries): " << babytot << endl;

	cout << "babywork: ";
	timestamp();

	mpf_set_z(dabs, discriminant);
	mpf_abs(dabs, dabs);
	mpf_sqrt(dsqrt, dabs);
	mpf_div(dsqrtoverpi, dsqrt, PIE);

	mpz_t leg_symbol, pmod;
	mpz_inits(leg_symbol, pmod, NULL);

	mpf_t Q, fone, fleg, accum, mainaccum, fp;
	mpf_inits(Q, fone, fleg, accum, mainaccum, fp, NULL);
	mpf_set(Q, dsqrtoverpi);
	mpf_set_ui(fone, 1);

	string line;
	// With prime36-38 is not delta representation.
	ifstream primefile("prime25.txt");
	bool DELTA_READ = true;

	assert(primefile.is_open());

	for (long i=0; i<6; i++) {
		mpf_init(ptasks[i].res);
	}

	// load in the prime table
	long p = 0;
	while (getline(primefile, line)) {

		if (DELTA_READ) {
			p += stol(line);
		} else {
			p = stol(line);
		}

		primetable[0] = p;

		long i = 1;
		for ( ; i < primetable.size() && getline(primefile, line); i++){
			if (DELTA_READ) {
				p += stol(line);
			} else {
				p = stol(line);
			}
			
			primetable[i] = p;
		}

		if (i != primetable.size()) {
			ptasks[0].lo = 0;
			ptasks[0].hi = i;
			primework(&ptasks[0]);
			mpf_mul(Q, Q, ptasks[0].res);
		} else {
			long prime_distance = primetable.size()/6;

			for (long j=0; j<6; j++) {
				ptasks[j].lo = prime_distance*j;
				if (j != 6-1) {
					ptasks[j].hi = prime_distance*(j+1);
				} else {
					ptasks[j].hi = primetable.size();
				}
			}

			for (int j=0; j<6; j++) {
				pthread_create(&slaves[j], NULL, primework, &ptasks[j]);
			}
			for (int j=0; j<6; j++) {
				pthread_join(slaves[j], NULL);
			}

			for(long j=0; j<6; j++) {
				mpf_mul(Q, Q, ptasks[j].res);
			}
		}
	}
	
	primefile.close();

	cout << "primework: ";
	timestamp();

	mpf_t margin;
	mpf_inits(margin, NULL);
	mpz_t margin_final, upper_bound, low_bound;
	mpz_inits(Q_final, margin_final, upper_bound, low_bound, NULL);

	mpz_set_f(Q_final, Q);

	mpf_sqrt(mainaccum, P);
	mpf_div(mainaccum, fone, mainaccum);
	mpf_mul(margin, Q, mainaccum);
	mpf_ceil(margin, margin);
	mpz_set_f(margin_final, margin);

	long HEURISTIC_CUT = 2; 
	mpz_tdiv_q_ui(margin_final, margin_final, HEURISTIC_CUT);

	mpz_add(upper_bound, Q_final, margin_final);
	mpz_sub(low_bound, Q_final, margin_final);

	mpf_t search_space_num, search_space_denom;
	mpf_inits(search_space_num, search_space_denom, NULL);
	mpf_set_z(search_space_num, margin_final);
	mpf_set_z(search_space_denom, Q_final);
	mpf_div(search_space_num, search_space_num, search_space_denom);
	cout << "paperbook search space ratio: " << search_space_num << endl;


	cout << "Q estimate: " << Q_final << endl;
	//cout << "lower bound estimate: " << low_bound << endl;
	//cout << "upper bound estimate: " << upper_bound << endl;

	//timestamp();

	mpz_init_set_ui(baby2, babytot*2);

	giant_step2 = form_copy(x);
	exponentiateinplace(giant_step2, baby2);

	form xp = form_copy(x);
	exponentiateinplace(xp, low_bound);

	mpz_t iter;
	mpz_init_set(iter, low_bound);

	mpz_t giantn;
	mpz_init(giantn);
	mpz_sub(giantn, upper_bound, low_bound);
	mpz_tdiv_q(giantn, giantn, baby2);

	// extend the upper bound so that it is the right termination condition for baby step giant step.
	mpz_add_ui(upper_bound, upper_bound, babytot);
	
	mpz_t task_distance, thread_steps;
	mpz_inits(task_distance, thread_steps, NULL);
	mpz_sub(task_distance, upper_bound, low_bound);
	mpz_tdiv_q_ui(task_distance, task_distance, 6);

	mpz_tdiv_q(thread_steps, task_distance, baby2);

	cout << "giant step no. iterations (per thread): " << thread_steps << endl;

	for (int i=0; i<6; i++) {

		mpz_init(gtasks[i].iter);
		mpz_mul_ui(gtasks[i].iter, task_distance, i);
		mpz_add(gtasks[i].iter, gtasks[i].iter, low_bound);

		mpz_init(gtasks[i].upper_bound);
		if (i != 6-1) {
			mpz_mul_ui(gtasks[i].upper_bound, task_distance, i+1); 
			mpz_add(gtasks[i].upper_bound, gtasks[i].upper_bound, low_bound);
		} else {
			mpz_set(gtasks[i].upper_bound, upper_bound);
		}
		
		mpz_inits(gtasks[i].xp.a, gtasks[i].xp.b, gtasks[i].xp.c, gtasks[i].xp.d, NULL);
		mpz_set(gtasks[i].xp.a, x.a);
		mpz_set(gtasks[i].xp.b, x.b);
		mpz_set(gtasks[i].xp.c, x.c);
		mpz_set(gtasks[i].xp.d, x.d);

		exponentiateinplace(gtasks[i].xp, gtasks[i].iter);

	}

	for (int i=0; i<6; i++) {
		pthread_create(&slaves[i], NULL, work, &gtasks[i]);
	}


	for (int i=0; i<6; i++) {
		pthread_join(slaves[i], NULL);
	}

	cout << "giantwork: ";
	timestamp();

	cout << "---- finished ----" << endl;

	return 0;
}
