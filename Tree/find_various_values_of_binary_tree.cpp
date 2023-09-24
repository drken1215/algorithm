//
// 強平衡二分木の各情報を求める
//   とくに、Frequency Table of Tree Distance を求める
//
// verified:
//   Codeforces Round 896 (Div. 1) C. Travel Plan
//     https://codeforces.com/contest/1868/problem/C
//
//   AtCoder ABC 321 E - Complete Binary Tree
//     https://atcoder.jp/contests/abc321/tasks/abc321_e
//


#include <bits/stdc++.h>
using namespace std;


// Find out of Strongly Balanced Binary Tree (N <= 10^18)
// the vertex number is 1-indexed (root = 1)
template<class mint> struct FindOutBinaryTree {
    // input data
    long long N;
    
    // main results
    vector<mint> depth_table;     // depth_table[d] := # of nodes whose distance from root is d
    vector<mint> distance_table;  // distance_tabls[l] := # of paths whose length is l
    
    // results of perfect binary trees
    vector<vector<mint>> perfect_depth_table, perfect_distance_table;
    
    // constructor
    FindOutBinaryTree() {}
    FindOutBinaryTree(long long n, bool build_dt = true) : N(n) {
        if (build_dt) init(n);
    }
    void set(long long n) {
        N = n;
    }
    void init(long long n) {
        N = n;
        int D = 0;
        while (n) { ++D, n /= 2; }
        findout_perfect_binary_tree(D);
        findout_binary_tree();
    }
    
    // preprocess of perfect binary trees
    void findout_perfect_binary_tree(int D) {
        auto pre = [&](auto self, long long d) -> vector<mint> {
            if (d == 0) {
                perfect_depth_table[d] = vector<mint>({mint(1)});
                return perfect_distance_table[d] = vector<mint>({mint(0), mint(1)});
            }
            vector<mint> depth(d+1, 0), distance(d*2+2, 0);
            for (int i = 0; i <= d; ++i) depth[i] = mint(1LL<<i);
            for (int i = 0; i <= d; ++i) distance[i+1] += mint(1LL<<i);
            for (int i = 1; i <= d; ++i) for (int j = 1; j <= d; ++j) {
                distance[i+j+1] += mint(1LL<<(i-1)) * mint(1LL<<(j-1));
            }
            const auto &left = self(self, d-1);
            for (int i = 0; i < left.size(); ++i) distance[i] += left[i] * 2;
            perfect_depth_table[d] = depth;
            return perfect_distance_table[d] = distance;
        };
        perfect_depth_table.resize(D+1);
        perfect_distance_table.resize(D+1);
        pre(pre, D);
    }
    
    // get left depth and right depth
    pair<long long, long long> get_depth(long long v) {
        long long left_depth = 0, right_depth = 0;
        long long left = v, right = v;
        while (left * 2 <= N) ++left_depth, left = left * 2;
        while (right * 2 + 1 <= N) ++right_depth, right = right * 2 + 1;
        return {left_depth, right_depth};
    }
    
    // find out the binary tree (size N)
    void findout_binary_tree() {
        auto rec = [&](auto self, long long v) -> pair<vector<mint>, vector<mint>> {
            vector<mint> depth, distance;
            if (v > N) return {depth, distance};
            
            // examine the depth of left subtree and right subtree
            auto [ld, rd] = get_depth(v);
            if (ld == rd) return {perfect_depth_table[ld], perfect_distance_table[rd]};

            // search the left subtree and right subtree
            auto [left_depth, left_distance] = self(self, v * 2);
            auto [right_depth, right_distance] = self(self, v * 2 + 1);
            depth.assign(max((int)left_depth.size(), (int)right_depth.size()) + 1, 0);
            distance.assign((int)left_depth.size() + (int)right_depth.size() + 2, 0);
            
            // update
            depth[0] = distance[1] = 1;
            for (int d = 0; d < (int)left_depth.size(); ++d) {
                depth[d + 1] += left_depth[d];
                distance[d + 2] += left_depth[d];
            }
            for (int d = 0; d < (int)right_depth.size(); ++d) {
                depth[d + 1] += right_depth[d];
                distance[d + 2] += right_depth[d];
            }
            for (int d1 = 0; d1 < (int)left_depth.size(); ++d1) {
                for (int d2 = 0; d2 < (int)right_depth.size(); ++d2) {
                    distance[d1 + d2 + 3] += left_depth[d1] * right_depth[d2];
                }
            }
            for (int l = 1; l < (int)left_distance.size(); ++l) {
                distance[l] += left_distance[l];
            }
            for (int l = 1; l < (int)right_distance.size(); ++l) {
                distance[l] += right_distance[l];
            }
            return {depth, distance};
        };
        auto [depth, distance] = rec(rec, 1);
        depth_table = depth;
        distance_table = distance;
    }
    
    // the number of nodes whose depth from v is d (v is 1-indexed)
    mint get_num_of_the_depth(long long v, long long d) {
        if (v <= 0 || v > N || d < 0) return mint(0);
        auto [left_depth, right_depth] = get_depth(v);
        if (left_depth < d) return mint(0);
        else if (right_depth >= d) return mint(1LL << d);
        else return mint(N - (v << d) + 1);
    }
    
    // the number of nodes whose distance from v is d (v is 1-indexed)
    mint get_num_of_the_distance(long long v, long long d) {
        if (v <= 0 || v > N) return mint(0);
        mint res = get_num_of_the_depth(v, d);
        for (long long i = 1; i <= d; ++i) {
            if (v == 1) break;
            if (i == d) {
                res += 1;
                break;
            }
            long long v2 = v / 2;
            if (v == v2 * 2 + 1) res += get_num_of_the_depth(v2 * 2, d - i - 1);
            else res += get_num_of_the_depth(v2 * 2 + 1, d - i - 1);
            v = v2;
        }
        return res;
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

void Codeforces_896_DIV1_C() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    auto solve = [&]() -> void {
        long long N, M;
        cin >> N >> M;
        FindOutBinaryTree<mint> fbt(N);
        mint res = 0;
        for (int l = 1; l < (int)fbt.distance_table.size(); ++l) {
            mint sum = 0;
            for (int i = 1; i < M; ++i) sum -= mint(i).pow(l);
            sum += mint(M).pow(l+1);
            sum *= mint(M).pow(N-l);
            res += fbt.distance_table[l] * sum;
        }
        cout << res << endl;
    };
    int T;
    cin >> T;
    while (T--) solve();
}

void ABC_321_E() {
    auto solve = [&]() -> void {
        long long N, X, D;
        cin >> N >> X >> D;
        FindOutBinaryTree<long long> fbt(N, false);
        cout << fbt.get_num_of_the_distance(X, D) << endl;
    };
    int T;
    cin >> T;
    while (T--) solve();
}


int main() {
    //Codeforces_896_DIV1_C();
    ABC_321_E();
}

