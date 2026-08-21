//
// Monge 最短路問題に関するアルゴリズム全集
//
// verified
//   AtCoder EDPC Z - Frog 3 (for simple LARSCH, LARSCH, monge CHT)
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry (for simple LARSCH, LARSCH, monge CHT)
//     https://codeforces.com/contest/319/problem/C 
// 
//   yukicoder No.705 ゴミ拾い Hard (for simple LARSCH, LARSCH, monge CHT)
//     https://yukicoder.me/problems/no/705 
//


#include <bits/stdc++.h>
using namespace std;


// LARSCH, in O(N)
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL, class FUNC> struct LARSCH {
    struct RowReduce;
    struct ColReduce;
    struct ColMap {
        const ColMap* parent = nullptr;
        const vector<int> *v = nullptr;
        constexpr int mapping(int i) const {
            int x = (v ? (*v)[i] : i);
            return (parent ? parent->mapping(x) : x);
        }
    };
    struct Eval {
        const FUNC *f = nullptr;
        long long a = 1, b = 0;  // row = a * i + b
        const ColMap *cm = nullptr;
        constexpr VAL operator () (int i, int j) const {
            int i2 = int(a * i + b), j2 = (cm ? cm -> mapping(j) : j);
            return (*f)(i2, j2);
        }
    };
    struct RowReduce {
        int N;
        Eval e;
        int cur_row = 0, state = 0;
        unique_ptr<ColReduce> rec;
        constexpr RowReduce(int N_, const Eval &e_) : N(N_), e(e_) {
            int M = N / 2;
            if (M > 0) {
                Eval e2 = e;
                e2.a = e.a * 2, e2.b = e.a + e.b;
                rec = make_unique<ColReduce>(M, e2);
            }
        }
        constexpr void reset() {
            cur_row = 0, state = 0;
            if (rec) rec->reset();
        }
        int get_argmin() {
            int i = cur_row++;
            if (!(i & 1)) {
                int prev = state, next = (i + 1 == N ? N - 1 : rec->get_argmin());
                state = next;
                int res = prev;
                for (int j = prev + 1; j <= next; j++) {
                    if (e(i, res) > e(i, j)) res = j;
                }
                return res;
            } else {
                return (e(i, state) <= e(i, i) ? state : i);
            }
        }
    };
    struct ColReduce {
        int N;
        Eval e;
        int cur_row = 0;
        vector<int> cols;
        ColMap cm;
        RowReduce rec;
        constexpr ColReduce(int N_, const Eval &e_) 
        : N(N_), e(e_), cols(), cm{e.cm, &cols}, rec(N_, Eval{e.f, e.a, e.b, &cm}) {
            cols.reserve(N);
        } 
        constexpr void reset() {
            cur_row = 0;
            cols.clear();
            rec.reset();
        }
        constexpr void push_col(int i, int j) {
            while (!cols.empty()) {
                int siz = (int)cols.size();
                if (siz == i) break;
                int last = cols.back();
                if (e(siz - 1, last) > e(siz - 1, j)) cols.pop_back();
                else break;
            }
            if ((int)cols.size() != N) cols.emplace_back(j);
        }
        constexpr int get_argmin() {
            int i = cur_row++;
            if (i == 0) {
                cols.clear();
                cols.emplace_back(0);
            } else {
                push_col(i, i * 2 - 1);
                push_col(i, i * 2);
            }
            return cols[rec.get_argmin()];
        }
    };

    FUNC f;
    ColMap root_cm;
    Eval root_eval;
    unique_ptr<RowReduce> base;
    explicit LARSCH(int N_, const FUNC &f_) 
    : f(std::move(f_)), root_cm{nullptr, nullptr}, root_eval{&f, 1, 0, &root_cm} {
        base = make_unique<RowReduce>(N_, root_eval);
    }
    constexpr void reset() {
        base->reset();
    }
    constexpr int get_argmin() {
        return base->get_argmin();
    }
};

// Monge Shortest Path Wrapper
template<class VAL> struct MongeShortestPath {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f) {
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;
        auto f2 = [&](int i, int j) -> VAL {
            i++;
            if (i <= j) return INF;
            else return dp[j] + f(j, i);
        };
        LARSCH<VAL, decltype(f2)> lar(N, f2);
        for (int r = 1; r <= N; r++) {
            int l = lar.get_argmin();
            dp[r] = dp[l] + f(l, r);
            prev[r] = l;
            cnt[r] = cnt[l] + 1;
        }
        res[0].first = VAL(0);
        for (int r = 1; r <= N; r++) res[r] = {dp[r], prev[r]};
        return res;
    }

    vector<int> reconstruct() {
        int N = (int)dp.size() - 1;
        vector<int> path;
        for (int v = N; v > 0; v = prev[v]) path.emplace_back(v);
        path.emplace_back(0);
        reverse(path.begin(), path.end());
        return path;
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    MongeShortestPath<long long> msp;
    const auto &res = msp.solve(N-1, func);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
/*
    A: 単調増加, B: 単調減少, ともに長さ N
    i -> j のコストが、B[i] × A[j] で与えられる ..... 単調増加 × 単調減少は Monge
    スタート: 0, ゴール: N-1
*/
void Codeforces_189_C() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    auto func = [&](int i, int j) -> long long {
        return B[i] * A[j];
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N-1, func);
    cout << res[N-1].first << endl;
}

// yukicoder No.705 ゴミ拾い Hard
/*
    A, X, Y: N 個
    これらを区間に分割していく
    　dp[j] = min_{0 ≦ i < j}(dp[i] + |A[j-1] - X[i]|^3 + |-Y[i]|^3)
    i -> j のコスト：|A[j-1] - X[i]|^3 + |-Y[i]|^3 ...... 差の凸関数 (Monge) + 縞々 (Monge) -> Monge
    スタート: 0, ゴール: N
*/
void yukicoder_705() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx * dx * dx + dy * dy * dy;
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N, func);
    cout << res[N].first << endl;
}


int main() {
    EDPC_Z();
    //Codeforces_189_C();
    //yukicoder_705();
}