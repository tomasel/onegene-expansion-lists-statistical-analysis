// From Prof. Blanzieri's lectures on Advanced Programming, University of Trento

#ifndef TESI_FRACTION_H
#define TESI_FRACTION_H

#include<iostream>
#include<set>
#include<cmath>

using namespace std;

class fraction {
    multiset<unsigned int> num, den;
    bool zero;
public:
    fraction() { zero = false; }

    fraction set_zero() {
        zero = true;
        num.clear();
        den.clear();
        return *this;
    }

    void mul(unsigned int _i);

    void div(unsigned int _i);

    fraction &operator*=(const fraction &_f);

    friend ostream &operator<<(ostream &os, const fraction &_f);

    long double compute2();     //versione forse più precisa

    bool is_zero() const { return zero; };
};

ostream &operator<<(ostream &os, const fraction &_f);

fraction hg(unsigned int _x, unsigned int _r, unsigned int _n, unsigned int _N);    // following wu's notation
fraction hg1(unsigned int x, unsigned int r, unsigned int n, unsigned int N);
fraction coeff_binom(unsigned int n, unsigned int k);                               //following wikipedia's notation

#endif //TESI_FRACTION_H
