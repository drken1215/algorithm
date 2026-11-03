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
//   ACL Contest 1 B - Sum is Multiple
//     https://atcoder.jp/contests/acl1/tasks/acl1_b
//


/*
    x ≡ b[i] (mod. m[i]) を満たす最小の正の整数を x として、x % MOD の値を求める
*/


#include <bits/stdc++.h>
using namespace std;


// Garner's algorithm
// for each step, we solve "coeffs[k] * t[k] + constants[k] = b[k] (mod. m[k])"
//      coeffs[k] = m[0]m[1]...m[k-1]
//      constants[k] = t[0] + t[1]m[0] + ... + t[k-1]m[0]m[1]...m[k-2]

// safe mod
template<class T_VAL, class T_MOD>
constexpr T_VAL safe_mod(T_VAL a, T_MOD m) {
    assert(m > 0);
    a %= m;
    if (a < 0) a += m;
    return a;
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

// if m is not coprime, call this function first
template<class T_VAL>
bool preGarner(vector<T_VAL> &b, vector<T_VAL> &m) {
    assert(b.size() == m.size());
    T_VAL res = 1;
    for (int i = 0; i < (int)b.size(); i++) {
        for (int j = 0; j < i; ++j) {
            T_VAL g = gcd(m[i], m[j]);
            if ((b[i] - b[j]) % g != 0) return false;
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
    vector<T_VAL> b2, m2;
    for (int i = 0; i < (int)b.size(); i++) {
        if (m[i] == 1) continue;
        b2.emplace_back(b[i]), m2.emplace_back(m[i]);
    }
    b = b2, m = m2;
    return true;
}

// find x, LCM (m must be coprime)
template<class T_VAL>
pair<T_VAL, T_VAL> Garner(const vector<T_VAL> &b, const vector<T_VAL> &m) {
    assert(b.size() == m.size());
    int num = (int)m.size();
    T_VAL res = 0, lcm = 1;
    vector<T_VAL> coeffs(num, 1), constants(num, 0);
    for (int k = 0; k < num; k++) {
        T_VAL t = safe_mod(b[k] - constants[k], m[k]) * mod_inv(coeffs[k], m[k]) % m[k];
        for (int i = k + 1; i < num; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
        res += t * lcm;
        lcm *= m[k];
    }
    return make_pair(res, lcm);
}

// find x (%MOD), LCM (%MOD) (m must be coprime)
template<class T_VAL, class T_MOD>
pair<T_VAL, T_VAL> Garner(const vector<T_VAL> &b, const vector<T_VAL> &m, T_MOD MOD) {
    assert(b.size() == m.size());
    assert(MOD > 0);
    int num = (int)m.size();
    T_VAL res = 0, lcm = 1;
    vector<T_VAL> coeffs(num, 1), constants(num, 0);
    for (int k = 0; k < num; k++) {
        T_VAL t = safe_mod(b[k] - constants[k], m[k]) * mod_inv(coeffs[k], m[k]) % m[k];
        for (int i = k + 1; i < num; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
        res = safe_mod(res + t * lcm, MOD);
        lcm = safe_mod(lcm * m[k], MOD);
    }
    return make_pair(res, lcm);
}


//------------------------------//
// Examples
//------------------------------//

// yukicoder 0187 中華風 (Hard)
void yukicoder_0187() {
    const int MOD = 1000000007;
    int N; 
    cin >> N;
    vector<long long> b(N), m(N);
    bool exist_non_zero = false;
    for (int i = 0; i < N; ++i) {
        cin >> b[i] >> m[i];
        if (b[i]) exist_non_zero = true;
    }
    if (!preGarner(b, m)) {
        cout << -1 << endl;
    } else if (!exist_non_zero) {
        long long lcm = 1;
        for (auto mi : m) lcm = lcm * mi % MOD;
        cout << lcm << endl;
    } else {
        auto [res, lcm] = Garner(b, m, MOD);
        cout << res << endl;
    }
}


// ACL Contest 1 B - Sum is Multiple
using i128 = __int128_t;
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
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
vector<long long> prime_factorize(long long n) {
    vector<long long> res;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        int num = 0;
        long long val = 1;
        while (n % p == 0) { ++num; n /= p; val *= p; }
        res.push_back(val);
    }
    if (n != 1) res.push_back(n);
    return res;
}
void ACL_Contest_1_B() {
    i128 N;
    cin >> N;
    N *= 2;
    auto pf = prime_factorize(N);
    int M = pf.size();
    i128 res = N - 1;
    for (int bit = 0; bit < (1<<M); ++bit) {
        i128 A = 1;
        for (int i = 0; i < M; ++i) if (bit & (1<<i)) A *= pf[i];
        i128 B = N / A;
        auto p = Garner(vector<i128>{0, B-1}, vector<i128>{A, B});
        if (p.first > 0) res = min(res, p.first);
        else res = min(res, p.second);
    }
    cout << res << endl;
}


int main() {
    yukicoder_0187();
    //ACL_Contest_1_B();
}