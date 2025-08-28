// Credits to Professor Blanzieri, from his lectures on Advanced Programming, University of Trento

#include "fraction.hpp"

using namespace std;

void fraction::mul(unsigned int _i) {
    auto it = den.find(_i);
    if (it != den.end()) den.erase(it);
    else num.insert(_i);
};


void fraction::div(unsigned int _i) {
    auto it = num.find(_i);
    if (it != num.end()) num.erase(it);
    else den.insert(_i);
}

fraction &fraction::operator*=(const fraction &_f) {
    if (_f.zero) {
        set_zero();
        return *this;
    }
    multiset<unsigned int>::iterator it;
    for (it = _f.num.begin(); it != _f.num.end(); it++) { mul(*it); }
    for (it = _f.den.begin(); it != _f.den.end(); it++) { div(*it); }
    return *this;
}

long double fraction::compute2() {
    if (zero) return 0;
    long double temp = 1;
    multiset<unsigned int>::reverse_iterator itn, itd;
    itn = num.rbegin();
    itd = den.rbegin();
    while ((itn != num.rend()) && (itd != den.rend())) {
		temp *= (double)*itn / (double)*itd;
        itd++;
        itn++;
    }
    while (itn != num.rend()) {
        temp *= *itn;
		if (itd != den.rend()) {
			temp /= *itd;
			++itd;
		}
        itn++;
    }
    while (itd != den.rend()) {
        temp /= *itd;
        itd++;
    }
    return temp;
}

ostream &operator<<(ostream &os, const fraction &_f) {
    if (_f.zero) return os << " 0/1 ";
    multiset<unsigned int>::iterator it;
    os << "( ";
    for (it = _f.num.begin(); it != _f.num.end(); it++) { os << *it << " "; }

    os << " ) / (";
    for (it = _f.den.begin(); it != _f.den.end(); it++) { os << " " << *it; }
    os << " )" << endl;
    return os;
}

fraction hg(unsigned int x, unsigned int r, unsigned int n, unsigned int N) {
    fraction temp;
    if ((((N > r + n) ? 0 : r - N + n) > x) || (((r < n) ? r : n) < x))
        return temp.set_zero();
    for (unsigned int i = n; i > 1; --i) temp.mul(i);
    for (unsigned int i = x; i > 1; --i) temp.div(i);
    for (unsigned int i = n - x; i > 1; --i) temp.div(i);
    for (unsigned int i = N - n; i > 1; --i) temp.mul(i);
    for (unsigned int i = r - x; i > 1; --i) temp.div(i);
    for (unsigned int i = N - n - r + x; i > 1; --i) temp.div(i);
    for (unsigned int i = N; i > 1; --i) temp.div(i);
    for (unsigned int i = r; i > 1; --i) temp.mul(i);
    for (unsigned int i = N - r; i > 1; --i) temp.mul(i);
    return temp;
}

fraction hg1(unsigned int x, unsigned int r, unsigned int n, unsigned int N) {
    fraction temp;
    if ((((N > r + n) ? 0 : r - N + n) > x) || (((r < n) ? r : n) < x))
        return temp.set_zero();
    for (unsigned int i = n; i > n - x; --i) temp.mul(i);
    for (unsigned int i = N - n; i > N - n - r + x; --i) temp.mul(i);
    for (unsigned int i = N; i > N - r; --i) temp.div(i);
    for (unsigned int i = x; i > 1; --i) temp.div(i);
    for (unsigned int i = r; i > r - x; --i) temp.mul(i);
    return temp;
}

fraction coeff_binom(unsigned int n, unsigned int k) {
    fraction temp;
    for (unsigned int i = n; i > n - k; --i) temp.mul(i);
    for (unsigned int i = k; i > 1; --i) temp.div(i);
    return temp;
}