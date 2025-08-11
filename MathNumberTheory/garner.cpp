//
// 中国剰余定理
//
// cf.
//   中国剰余定理 (CRT) の解説と、それを用いる問題のまとめ
//     https://qiita.com/drken/items/ae02240cd1f8edfc86fd
//
// verified
//   yukicoder 0187 中華風 (Hard)
//     https://yukicoder.me/problems/no/187
//


/*
    x ≡ b[i] (mod. m[i]) を満たす最小の正の整数を x として、x % MOD の値を求める
*/


#include <bits/stdc++.h>
using namespace std;


// safe mod
template<class T_VAL, class T_MOD>
constexpr T_VAL safe_mod(T_VAL a, T_MOD m) {
    assert(m > 0);
    a %= m;
    if (a < 0) a += m;
    return a;
}

// mod pow
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_pow(T_VAL a, T_VAL n, T_MOD m) {
    assert(m > 0);
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
    assert(m > 0);
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

// Garner's algorithm
// if m is not coprime, call this function first
template<class T_VAL, class T_MOD>
T_VAL preGarner(vector<T_VAL> &b, vector<T_VAL> &m, T_MOD MOD) {
    T_VAL res = 1;
    for (int i = 0; i < (int)b.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            T_VAL g = gcd(m[i], m[j]);
            if ((b[i] - b[j]) % g != 0) return -1;
            m[i] /= g, m[j] /= g;
            T_VAL gi = gcd(m[i], g), gj = g/gi;
            do {
                g = gcd(gi, gj);
                gi *= g, gj /= g;
            } while (g != 1);
            m[i] *= gi, m[j] *= gj;
            b[i] %= m[i], b[j] %= m[j];
        }
    }
    for (int i = 0; i < (int)b.size(); ++i) (res *= m[i]) %= MOD;
    return res;
}

// x (%MOD), LCM (%MOD) を求める (m must be coprime)
// for each step, we solve "coeffs[k] * t[k] + constants[k] = b[k] (mod. m[k])"
//      coeffs[k] = m[0]m[1]...m[k-1]
//      constants[k] = t[0] + t[1]m[0] + ... + t[k-1]m[0]m[1]...m[k-2]
template<class T_VAL, class T_MOD>
T_VAL Garner(vector<T_VAL> b, vector<T_VAL> m) {
    assert(b.size() == m.size());
    int num = (int)m.size();
    vector<T_VAL> coeffs(num + 1, T_VAL(1)), constants(num + 1, T_VAL(0));
    for (int k = 0; k < num; k++) {
        T_VAL t = safe_mod(safe_mod(b[k] - constants[k], m[k]) * mod_inv(coeffs[k], m[k]), m[k]);
        for (int i = k + 1; i < num; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
        constants.back() += t * coeffs.back();
        coeffs.back() *= m[k];
    }
    return constants.back();
}

template<class T_VAL, class T_MOD>
T_VAL Garner(vector<T_VAL> b, vector<T_VAL> m, T_MOD MOD) {
    assert(b.size() == m.size());
    assert(MOD > 0);
    int num = (int)m.size();
    vector<T_VAL> coeffs(num + 1, T_VAL(1)), constants(num + 1, T_VAL(0));
    m.emplace_back(MOD);  // banpei
    for (int k = 0; k < num; k++) {
        T_VAL t = safe_mod(safe_mod(b[k] - constants[k], m[k]) * mod_inv(coeffs[k], m[k]), m[k]);
        for (int i = k + 1; i < num + 1; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
    }
    return constants.back();
}


//------------------------------//
// Examples
//------------------------------//

// yukicoder 0187 中華風 (Hard)
void yukicoder_0187() {
    int N; 
    cin >> N;
    vector<long long> b(N), m(N);
    bool exist_non_zero = false;
    for (int i = 0; i < N; ++i) {
        cin >> b[i] >> m[i];
        if (b[i]) exist_non_zero = true;
    }
    const int MOD = 1000000007;
    long long lcm = preGarner(b, m, MOD);
    
    if (!exist_non_zero) cout << lcm << endl;
    else if (lcm == -1) cout << -1 << endl;
    else cout << Garner(b, m, MOD) << endl;
}


int main() {
    yukicoder_0187();
}