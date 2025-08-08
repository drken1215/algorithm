//
// mod_pow, mod_inv
//
// reference:
//   drken: 「998244353 で割ったあまり」の求め方を総特集！ ～ 逆元から離散対数まで ～
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   simple test
//


#include <bits/stdc++.h>
using namespace std;


// mod pow
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_pow(T_VAL a, T_VAL n, T_MOD m) {
    T_VAL res = 1;
    while (n > 0) {
        if (n % 2 == 1) res = res * a % m;
        a = a * a % m;
        n >>= 1;
    }
    return res;
}

// mod inv
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_inv(T_VAL a, T_MOD m) {
    T_VAL b = m, u = 1, v = 0;
    while (b > 0) {
        T_VAL t = a / b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= m;
    if (u < 0) u += m;
    return u;
}



//------------------------------//
// Examples
//------------------------------//

using i128 = __int128_t;
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
         --d;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        --d;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}

void small_check() {

    const i128 P = 9999999999971LL;
    for (i128 a = 1; a <= 100; ++a) {
        cout << mod_pow(a, P - 2, P) << ", " << mod_inv(a, P) << endl;
    }
}


int main() {
    small_check();
}