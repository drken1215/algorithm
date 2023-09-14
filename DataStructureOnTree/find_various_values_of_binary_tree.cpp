//
// 強平衡二分木の各情報を求める
//
// verified:
//   Codeforces Round 896 (Div. 1) C. Travel Plan
//     https://codeforces.com/contest/1868/problem/C
//


#include <bits/stdc++.h>
using namespace std;


// Find out of Strongly Balanced Binary Tree (size <= 10^18)
template<class mint> struct FindOutBinaryTree {
    // input data
    long long N;
    
    // main results
    vector<mint> num_nodes;  // num_nodes[d] := # of nodes whose distance from root is d
    vector<mint> num_paths;  // num_paths[l] := # of paths whose length is l
    
    // results of perfect binary trees
    vector<vector<mint>> perfect_num_nodes, perfect_num_paths;
    
    // constructor
    FindOutBinaryTree() {}
    FindOutBinaryTree(long long n) : N(n) { init(n); }
    void init(long long n) {
        N = n;
        int D = 0;
        while (n) { ++D, n /= 2; }
        preprocess_perfect_binary_tree(D);
        findout_binary_tree();
    }
    
    // preprocess of perfect binary trees
    void preprocess_perfect_binary_tree(int D) {
        auto pre = [&](auto self, long long d) -> vector<mint> {
            if (d == 0) {
                perfect_num_nodes[d] = vector<mint>({mint(1)});
                return perfect_num_paths[d] = vector<mint>({mint(0), mint(1)});
            }
            vector<mint> num_nodes(d+1, 0), num_paths(d*2+2, 0);
            for (int i = 0; i <= d; ++i) num_nodes[i] = mint(1LL<<i);
            for (int i = 0; i <= d; ++i) num_paths[i+1] += mint(1LL<<i);
            for (int i = 1; i <= d; ++i) for (int j = 1; j <= d; ++j) {
                num_paths[i+j+1] += mint(1LL<<(i-1)) * mint(1LL<<(j-1));
            }
            const auto &left = self(self, d-1);
            for (int i = 0; i < left.size(); ++i) num_paths[i] += left[i] * 2;
            perfect_num_nodes[d] = num_nodes;
            return perfect_num_paths[d] = num_paths;
        };
        perfect_num_nodes.resize(D+1);
        perfect_num_paths.resize(D+1);
        pre(pre, D);
    }
    
    // find out the binary tree (size N)
    void findout_binary_tree() {
        auto rec = [&](auto self, long long v) -> pair<vector<mint>, vector<mint>> {
            vector<mint> num, dp;
            if (v > N) return {num, dp};
            
            int left_depth = 0, right_depth = 0;
            long long left_val = v, right_val = v;
            while (left_val * 2 <= N) ++left_depth, left_val = left_val * 2;
            while (right_val * 2 + 1 <= N) ++right_depth, right_val = right_val * 2 + 2;
                        
            // 完全二分木の場合
            if (left_depth == right_depth) {
                return {perfect_num_nodes[left_depth], perfect_num_paths[left_depth]};
            }

            // 左右の二分木を探索
            auto [left_num, left_dp] = self(self, v * 2);
            auto [right_num, right_dp] = self(self, v * 2 + 1);
            num.assign(left_depth + 1, 0);
            dp.assign((int)left_num.size() + (int)right_num.size() + 2, 0);
            
            // 更新
            num[0] = 1;
            dp[1] = 1;
            for (int d = 0; d < (int)left_num.size(); ++d) {
                num[d + 1] += left_num[d];
                dp[d + 2] += left_num[d];
            }
            for (int d = 0; d < (int)right_num.size(); ++d) {
                num[d + 1] += right_num[d];
                dp[d + 2] += right_num[d];
            }
            for (int d1 = 0; d1 < (int)left_num.size(); ++d1) {
                for (int d2 = 0; d2 < (int)right_num.size(); ++d2) {
                    dp[d1 + d2 + 3] += left_num[d1] * right_num[d2];
                }
            }
            for (int l = 1; l < (int)left_dp.size(); ++l) dp[l] += left_dp[l];
            for (int l = 1; l < (int)right_dp.size(); ++l) dp[l] += right_dp[l];
            return {num, dp};
        };
        auto [num, dp] = rec(rec, 1);
        num_nodes = num;
        num_paths = dp;
    }
};

// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() noexcept : val(0) { }
    constexpr Fp(long long v) noexcept : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr long long get() const noexcept { return val; }
    constexpr int get_mod() const noexcept { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator - () const noexcept {
        return val ? MOD - val : 0;
    }
    constexpr Fp operator + (const Fp &r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) noexcept {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) noexcept {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) noexcept {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) noexcept {
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
    constexpr Fp pow(long long n) const noexcept {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const noexcept {
        Fp res(1), div(*this);
        return res / div;
    }

    // other operators
    constexpr bool operator == (const Fp &r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) noexcept {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD> &r, long long n) noexcept {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> modinv(const Fp<MOD> &r) noexcept {
        return r.inv();
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void Codeforces_896_DIV1_C_Solve() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    long long N, M;
    cin >> N >> M;
    
    FindOutBinaryTree<mint> fbt(N);
    mint res = 0;
    for (int l = 1; l < (int)fbt.num_paths.size(); ++l) {
        mint sum = 0;
        for (int i = 1; i < M; ++i) sum -= mint(i).pow(l);
        sum += mint(M).pow(l+1);
        sum *= mint(M).pow(N-l);
        
        res += fbt.num_paths[l] * sum;
    }
    cout << res << endl;
}

void Codeforces_896_DIV1_C() {
    int T;
    cin >> T;
    while (T--) Codeforces_896_DIV1_C_Solve();
}


int main() {
    Codeforces_896_DIV1_C();
}


