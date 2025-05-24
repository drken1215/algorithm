//
// RMQ (by sparse table)
//
// verified:
//   ABC 407 F - Sums of Sliding Window Maximum
//     https://atcoder.jp/contests/abc407/tasks/abc407_f
//


/*
  動的更新はできないが、最初に静的に構築してしまえば、範囲最小値を O(1) でとり出せる
    交叉半束 (meet-semi-lattice) 上で動作
    半群 (semi-group) 上で動作する上位互換バージョン (disjoint sparse table) もある

  init(vec): 配列を vec で初期化構築, O(n logn)
  get(a, b): [a, b) の最小値を取得
*/


#include <bits/stdc++.h>
using namespace std;


// Sparse Table
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP;
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() { }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  { init(vec, op); }
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<MeetSemiLattice>(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            for (int j = 0; j < n; ++j)
                dat[i][j] = OP(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
        }
    }
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 407 F - Sums of Sliding Window Maximum
void ABC_407_F() {
    int N;
    cin >> N;
    using Node = pair<long long, int>;
    vector<Node> A(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i].first;
        A[i].second = i;
    }
    SparseTable<Node> st(A, [&](Node a, Node b){return max(a, b);});

    // imos method
    vector<long long> res(N+3, 0);
    auto imos = [&](int l, int r, long long fval, long long d) -> void {
        // 区間 [l, r) に、初項 fval, 公差 d の等差数列を足す
        long long lval = fval + d * (r - l - 1);
        res[l] += fval, res[l+1] += -fval + d, res[r] += -lval - d, res[r+1] += lval;
    };

    // 再帰的に最大値を足し込む
    auto rec = [&](auto &&rec, int l, int r) -> void {
        if (r - l <= 0) return;
        auto [val, m] = st.get(l, r);
        int L = min(m+1-l, r-m), R = max(m+1-l, r-m);
        imos(1, L+1, val, val);
        imos(L+1, R+1, val * L, 0);
        imos(R+1, r-l+1, val * (L-1), -val);
        rec(rec, l, m);
        rec(rec, m+1, r);
    };
    rec(rec, 0, N);

    // imos の復元 (2 回累積和をとる)
    for (int i = 0; i+1 < res.size(); i++) res[i+1] += res[i];
    for (int i = 0; i+1 < res.size(); i++) res[i+1] += res[i];

    // 出力
    for (int k = 1; k <= N; k++) cout << res[k] << '\n';
}


int main() {
    ABC_407_F();
}