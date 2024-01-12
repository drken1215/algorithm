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


#include <iostream>
#include <vector>
#include <bitset>
using namespace std;


const long long MOD = 1000000007;

// 最大公約数
long long GCD(long long a, long long b) {
    if (b == 0) return a;
    else return GCD(b, a % b);
}

// Garner のアルゴリズムの前処理
long long PreGarner(vector<long long> &b, vector<long long> &m, long long MOD) {
    long long res = 1;
    for (int i = 0; i < (int)b.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            long long g = GCD(m[i], m[j]);
            if ((b[i] - b[j]) % g != 0) return -1;
            m[i] /= g;
            m[j] /= g;
            long long gi = GCD(m[i], g);
            long long gj = g/gi;
            do {
                g = GCD(gi, gj);
                gi *= g, gj /= g;
            } while (g != 1);
            m[i] *= gi, m[j] *= gj;
            b[i] %= m[i], b[j] %= m[j];
        }
    }
    for (int i = 0; i < (int)b.size(); ++i) (res *= m[i]) %= MOD;
    return res;
}

// 負の数にも対応した mod (a = -11 とかでも OK)
inline long long mod(long long a, long long m) {
    long long res = a % m;
    if (res < 0) res += m;
    return res;
}

// 拡張 Euclid の互除法
long long extGCD(long long a, long long b, long long &p, long long &q) {
    if (b == 0) { p = 1; q = 0; return a; }
    long long d = extGCD(b, a%b, q, p);
    q -= a/b * p;
    return d;
}

// 逆元計算 (ここでは a と m が互いに素であることが必要)
long long modinv(long long a, long long m) {
    long long x, y;
    extGCD(a, m, x, y);
    return mod(x, m); // 気持ち的には x % m だが、x が負かもしれないので
}

// Garner のアルゴリズム, x%MOD, LCM%MOD を求める (m は互いに素でなければならない)
// for each step, we solve "coeffs[k] * t[k] + constants[k] = b[k] (mod. m[k])"
//      coeffs[k] = m[0]m[1]...m[k-1]
//      constants[k] = t[0] + t[1]m[0] + ... + t[k-1]m[0]m[1]...m[k-2]
long long Garner(vector<long long> b, vector<long long> m, long long MOD) {
    m.push_back(MOD); // banpei
    vector<long long> coeffs((int)m.size(), 1);
    vector<long long> constants((int)m.size(), 0);
    for (int k = 0; k < (int)b.size(); ++k) {
        long long t = mod((b[k] - constants[k]) * modinv(coeffs[k], m[k]), m[k]);
        for (int i = k+1; i < (int)m.size(); ++i) {
            (constants[i] += t * coeffs[i]) %= m[i];
            (coeffs[i] *= m[k]) %= m[i];
        }
    }
    return constants.back();
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int N; cin >> N;
    vector<long long> b(N), m(N);
    bool exist_non_zero = false;
    for (int i = 0; i < N; ++i) {
        cin >> b[i] >> m[i];
        if (b[i]) exist_non_zero = true;
    }
    long long lcm = PreGarner(b, m, MOD);
    
    if (!exist_non_zero) cout << lcm << endl;
    else if (lcm == -1) cout << -1 << endl;
    else cout << Garner(b, m, MOD) << endl;
}
