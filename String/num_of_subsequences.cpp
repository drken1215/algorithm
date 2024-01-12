//
// 文字列 S の相異なる部分列の個数を求める
//
// reference:
//   drken: 部分列 DP --- 文字列の部分文字列を重複なく走査する DP の特集
//     https://qiita.com/drken/items/a207e5ae3ea2cf17f4bd
//
// verified:
//   AtCoder TDPC  G - 辞書順
//     https://atcoder.jp/contests/tdpc/tasks/tdpc_lexicographical
//


#include <bits/stdc++.h>
using namespace std;


// dp[i][c] := #subsequences of S.substr(i) whose first letter is c
// nex[i][c] := the first index of letter c in S.substr(i)
template<class T> struct RunSubsequences {
    // inner value
    int N;
    string S;
    vector<vector<T>> dp;
    vector<vector<int>> nex;
    
    // constructors
    RunSubsequences() {}
    RunSubsequences(const string &s, bool calc_nex = true) : S(s) {
        init(s, calc_nex);
    }
    
    // initializer
    void init(const string &s, bool calc_nex = true) {
        S = s;
        N = (int)S.size();
        init_dp();
        if (calc_nex) init_nex();
    }
    void init_dp() {
        dp.assign(N+1, vector<T>(26, 0));
        dp[N-1][S[N-1]-'a'] = 1;
        for (int i = N-2; i >= 0; --i) {
            for (int c = 0; c < 26; ++c) {
                if (c != S[i]-'a') dp[i][c] += dp[i+1][c];
                else {
                    dp[i][c] += 1;  // only "c"
                    for (int c2 = 0; c2 < 26; ++c2) dp[i][c] += dp[i+1][c2];
                }
            }
        }
    }
    void init_nex() {
        nex.assign(N+2, vector<int>(26, N+1));
        for (int i = N-1; i >= 0; --i) {
            for (int j = 0; j < 26; ++j) nex[i][j] = nex[i+1][j];
            nex[i][S[i]-'a'] = i;
        }
    }
    
    // k-th smallest subsequence
    string reconstruct(T K) {
        T sum = 0;
        for (int c = 0; c < 26; ++c) sum += dp[0][c];
        if (sum < K) return "";
        
        string res = "";
        for (int i = 0; i < N; ++i) {
            char c;
            for (int j = 0; j < 26; ++j) {
                if (K - dp[i][j] <= 0) {
                    c = 'a'+j;
                    break;
                }
                K -= dp[i][j];
            }
            res += c;
            K -= 1;
            if (K <= 0) break;
            while (S[i] != c) ++i;
        }
        return res;
    }
};


// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() : val(0) { }
    constexpr Fp(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    
    // getter
    constexpr long long get() const {
        return val;
    }
    constexpr int get_mod() const {
        return MOD;
    }
    
    // comparison operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        val = val * u % MOD;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp(0) - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    
    // other operators
    constexpr Fp& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) {
        return os << x.val;
    }
    
    // other functions
    constexpr Fp pow(long long n) const {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        Fp res(1), div(*this);
        return res / div;
    }
    friend constexpr Fp<MOD> pow(const Fp<MOD> &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> inv(const Fp<MOD> &r) {
        return r.inv();
    }
};



//------------------------------//
// Examples
//------------------------------//

void TDPC_G() {
    // 上限が INF になる long long 型
    const long long INF = 1LL<<60;
    struct Val {
        long long val;
        Val(long long v = 0) : val(v) {}
        Val& operator += (const Val &r) {
            val += r.val;
            if (val >= INF) val = INF;
            return *this;
        }
        Val& operator -= (const Val &r) {
            val -= r.val;
            return *this;
        }
        Val operator + (const Val &r) const { return Val(*this) += r; }
        Val operator - (const Val &r) const { return Val(*this) -= r; }
        bool operator < (const Val &r) const { return val < r.val; }
        bool operator <= (const Val &r) const { return val <= r.val; }
    };
    
    string S;
    long long K;
    cin >> S >> K;
    RunSubsequences<Val> rs(S, false);
    string res = rs.reconstruct(K);
    cout << (res != "" ? res : "Eel") << endl;
}


int main() {
    TDPC_G();
}
