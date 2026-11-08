
//
// Monge 最短路問題に関するアルゴリズム全集
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君 (for Monge CHT)
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3 (for simple LARSCH, LARSCH, monge CHT)
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry (for simple LARSCH, LARSCH, monge CHT)
//     https://codeforces.com/contest/319/problem/C 
// 
//   yukicoder No.705 ゴミ拾い Hard (for simple LARSCH, LARSCH, monge CHT)
//     https://yukicoder.me/problems/no/705 
//
//   ABC 218 H - Red and Blue Lamps (for Alien DP)
//     https://atcoder.jp/contests/abc218/tasks/abc218_h
//
//   yukicoder No.952 危険な火薬庫 (for SMAWK + D-enum)
//     https://yukicoder.me/problems/no/952
//
//   Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
//     2-solutions (slide access Monotone Minima + D-enum, slide access simple LARSCH + Alien DP)
//     https://codeforces.com/contest/868/problem/F
//
//   ABC 480 G - Zonal Score Minimization (for anti-Monge D&D SMAWK, Monge CHT)
//     https://atcoder.jp/contests/abc480/tasks/abc480_g
//


#include <bits/stdc++.h>
using namespace std;



/*
    Alien's Trick (by ラグランジュ緩和)

    min_{x}: f(x) s.t. g(x) = K
    = max_{λ は整数}: min_{x} (f(x) + λ(g(x) - K))

    max_{x}: f(x) s.t. g(x) = K
    = min_{λ は整数}: max_{x} (f(x) + λ(g(x) - K))

    【仮定】
　　・f(x), g(x): 0 以上の整数値をとる関数
　　・K: 整数
　　・g(x) = p となる x についての f(x) の最小値を h(p) とおくと、h(p) は p について下に凸な関数
　　・h(K) < ∞ (g(x) = K となる x が存在)

    【インターフェース】
　　・lag(λ): x✳︎ = argmin_{x} {f(x) + λ(g(x) - K)} としたときのペア値 {f(x✳︎), g(x✳︎)} を返す関数
　　・返り値: {最適値, そのときの λ}
*/
template<class VAL> struct AliensTrick {
    // find lo such that g(lo) >= R
    template<class FUNC> tuple<VAL, VAL, VAL> find_low(const FUNC &lag, VAL R) {
        auto [f0, g0] = lag(VAL(0));
        VAL lo = 0, hi = 0, flo = f0, glo = g0;
        if (glo >= R) return tie(lo, flo, glo);
        lo = VAL(-1);
        tie(flo, glo) = lag(lo);
        while (glo < R) {
            assert(lo < hi);
            VAL dif = hi - lo;
            hi = lo, lo -= VAL(2) * dif;
            tie(flo, glo) = lag(lo);
        }
        return tie(lo, flo, glo);
    }

    // find hi such that g(hi) <= L
    template<class FUNC> tuple<VAL, VAL, VAL> find_high(const FUNC &lag, VAL L) {
        auto [f0, g0] = lag(VAL(0));
        VAL lo = 0, hi = 0, fhi = f0, ghi = g0;
        if (ghi <= L) return tie(hi, fhi, ghi);
        hi = VAL(1);
        tie(fhi, ghi) = lag(hi);
        while (ghi > L) {
            assert(lo < hi);
            VAL dif = hi - lo;
            lo = hi, hi += VAL(2) * dif;
            tie(fhi, ghi) = lag(hi);
        }
        return tie(hi, fhi, ghi);
    }

    // solver
    template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL K) {
        // zero check
        auto [f0, g0] = lag(VAL(0));
        if (g0 == K) return {f0, VAL(0)};

        // find an initial value of binary search (by doubling)
        auto [lo, flo, glo] = find_low(lag, K);
        auto [hi, fhi, ghi] = find_high(lag, K);
        if (glo == K) return {flo, lo};
        if (ghi == K) return {fhi, hi};
        assert(glo > K && ghi < K);

        // binary search
        while (hi - lo > VAL(1)) {
            VAL mid = (lo + hi) / 2;
            auto [fmid, gmid] = lag(mid);
            if (gmid == K) return {fmid, mid};
            if (gmid > K) lo = mid, flo = fmid, glo = gmid;
            else hi = mid, fhi = fmid, ghi = gmid;
        }

        // in case: h(p) has a collinear part (= could not achieve g(x) = K)
        // return the Lagrange value
        VAL reslo = flo + lo * (glo - K), reshi = fhi + hi * (ghi - K);
        if (reslo >= reshi) return {reslo, lo};
        else return {reshi, hi};
    }
    template<class FUNC> pair<VAL, VAL> minimize(const FUNC &lag, VAL K) {
        return solve(lag, K);
    }
    template<class FUNC> pair<VAL, VAL> maximize(const FUNC &lag, VAL K) {
        auto reverse_lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto [f, g] = lag(-lambda);
            return {-f, g};
        };
        auto [f, g] = solve(reverse_lag, K);
        return {-f, -g};
    }

    // in case: L <= g(x) <= R
    template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL L, VAL R) {
        assert(L <= R);
        auto [f0, g0] = lag(VAL(0));
        if (L <= g0 && g0 <= R) return {f0, VAL(0)};
        if (g0 > R) return solve(lag, R);
        else return solve(lag, L);
    }
    template<class FUNC> pair<VAL, VAL> minimize(const FUNC &lag, VAL L, VAL R) {
        return solve(lag, L, R);
    }
    template<class FUNC> pair<VAL, VAL> maximize(const FUNC &lag, VAL L, VAL R) {
        auto reverse_lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto [f, g] = lag(-lambda);
            return {-f, g};
        };
        auto [f, g] = solve(reverse_lag, L, R);
        return {-f, -g};
    }
};

// Alien DP
// find the d-edges shortest path from vertex 0 on DAG by Lagrange relaxation
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct AlienDP {
    // inner values
    MongeShortestPath<VAL> msp;
    AliensTrick<VAL> at;

    // solver (random access ver.)
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, int D
    , const string &solver = "simple_larsch") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &obj = msp.solve(N, fg, solver);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

    // solver (slide access ver.)
    template<class STATE, class ADD, class DEL, class GETCOST> pair<VAL, int> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, int D
    , const string &solver = "simple_larsch") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto get_lag = [&](STATE &st) -> VAL { return get(st) + lambda; };
            const auto &obj = msp.solve(N, ini, add, del, get_lag, solver);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

    // debugger
    template<class FUNC> void debug(int N, const FUNC &f, VAL start = 0, VAL goal = 10) {
        auto lag = [&](VAL lambda) -> void {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &dp = msp.solve(N, fg);
            auto cnt = msp.cnt[N];
            auto res = dp.back().first - lambda * cnt;
            cout << lambda << ": " << res << "(" << cnt << "); (";
            for (int i = 0; i < (int)dp.size(); i++) {
                if (i) cout << " ";
                cout << "<" << dp[i].first << ", " << dp[i].second << ">";
            }
            cout << ")" << endl;
        };
        for (VAL lambda = start; lambda <= goal; lambda++) lag(lambda);
    }
};

// f: Monge, interval monotone ((N+1) x (N+1))
// dp[i][j] = min_{i <= k < j} (dp[i][k] + dp[k][j]) + f(i, j), in O(N^2)
template<class VAL, class FUNC> VAL KnuthYao(int N, const FUNC &f) {
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<vector<VAL>> dp(N+1, vector<VAL>(N+1, INF));
    vector<vector<int>> K(N+1, vector<int>(N+1, -1));

    for (int i = 0; i <= N; i++) {
        dp[i][i] = 0, K[i][i] = i;
        if (i < N) dp[i][i+1] = f(i, i+1), K[i][i+1] = i;
    }
    for (int between = 2; between <= N; between++) {
        for (int i = 0; i + between <= N; i++) {
            int j = i + between;
            for (int k = K[i][j-1]; k <= K[i+1][j]; k++) {
                VAL tmp = dp[i][k] + dp[k][j] + f(i, j);
                if (dp[i][j] > tmp) {
                    dp[i][j] = tmp;
                    K[i][j] = k;
                }   
            }
        }
    }
    return dp[0][N];
}


//------------------------------//
// Examples
//------------------------------//

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    DynamicLiChaoTreeForPSC<long long> cht(0, 210000, 1LL<<60);
    for (long long i = 0; i < N; ++i) {
        auto func = [i, &a](long long x) -> long long { return (x - i) * (x - i) + a[i]; };
        cht.insert(func);
    }
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i).first;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
void EDPC_Z(const string solver = "simple_larsch") {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    MongeShortestPath<long long> msp;
    const auto &res = msp.solve(N-1, func, solver);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
/*
    A: 単調増加, B: 単調減少, ともに長さ N
    i -> j のコストが、B[i] × A[j] で与えられる ..... 単調増加 × 単調減少は Monge
    スタート: 0, ゴール: N-1
*/
void Codeforces_189_C(const string solver = "simple_larsch") {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    auto func = [&](int i, int j) -> long long {
        return B[i] * A[j];
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N-1, func, solver);
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
void yukicoder_705(const string solver = "simple_larsch") {
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
    auto res = msp.solve(N, func, solver);
    cout << res[N].first << endl;
}

// ABC 218 H - Red and Blue Lamps
void ABC_218_H() {
    long long N, R;
    cin >> N >> R;
    vector<long long> A(N - 1);
    for (int i = 0; i < N - 1; i++) cin >> A[i];
    auto cost = [&](int i, int j) -> long long {
        if (j - i <= 1) return 0;
        if (i == 0) return -A[j - 2];
        else return -A[i - 1] - A[j - 2];
    };
    AlienDP<long long> ad;
    //ad.debug(N, cost);
    auto [resR, cntR] = ad.solve(N, cost, R);  // both "larsch" and "simple_larsch" are OK
    auto [resB, cntB] = ad.solve(N, cost, N - R);
    auto res = max(-resR, -resB);
    cout << res << endl;
}

// yukicoder No.952 危険な火薬庫
/*
    末尾にドアを追加する。
    ノード 0, 1, ..., N, N+1
    f(i, j) = (S[j-1] - S[i])^2
    選ぶドアが k 個 → 区間の個数は、N+1-k
*/
void yukicoder_952() {
    int N;
    cin >> N;
    vector<long long> A(N + 1, 0), S(N + 2, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i + 1] = S[i] + A[i];
    auto f = [&](int i, int j) -> long long { 
        return (S[j - 1] - S[i]) * (S[j - 1] - S[i]);
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N + 1, f, N + 1);
    for (int k = 1; k <= N; k++) cout << res[N + 1 - k] << endl;
}

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Monge 単一始点 d-辺最短路の d = 1, 2, ..., D における列挙 (by SMAWK 法, in O(ND))
void Codeforces438_F_enum() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    vector<long long> nums(N + 1, 0);
    long long left = 0, right = 0, score = 0;
    auto add = [&](int i) -> void {
        score += nums[A[i]];
        nums[A[i]]++;
    };
    auto del = [&](int i) -> void {
        nums[A[i]]--;
        score -= nums[A[i]];
    };
    auto f = [&](int i, int j) -> long long {
        while (left < i) del(left++);
        while (i < left) add(--left);
        while (right < j) add(right++);
        while (j < right) del(--right);
        return score;
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N, f, K, "monotone");  // we should use slide access
    cout << res[K] << endl;
}

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Alien DP (not the assumed solution)
void Codeforces438_F_alien() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    struct STATE {
        vector<long long> nums;
        long long cost;
        STATE(int N) : nums(N + 1, 0), cost(0) {}
    };
    STATE ini(N);
    auto add = [&](STATE &st, int i) -> void {
        st.cost += st.nums[A[i]];
        st.nums[A[i]]++;
    };
    auto del = [&](STATE &st, int i) -> void {
        st.nums[A[i]]--;
        st.cost -= st.nums[A[i]];
    };
    auto get = [&](STATE &st) -> long long {
        return st.cost;
    };

    AlienDP<long long> adp;
    auto res = adp.solve(N, ini, add, del, get, K);  // we should use slide access
    cout << res.first << endl;
}

// ABC 480 G - Zonal Score Minimization
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP = [](const MeetSemiLattice &l, const MeetSemiLattice &r) {
        return min(l, r);
    };
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() {}
    SparseTable(const vector<MeetSemiLattice> &vec) {
        init(vec);
    }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  {
        init(vec, op);
    }
    void init(const vector<MeetSemiLattice> &vec) {
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
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
        init(vec);
    }
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};
void ABC_480_G_SMAWK() {
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});

    auto func = [&](int l, int r) -> long long {
        assert(r - l >= 2);
        return stma.get(l, r) - stmi.get(l, r);
    };
    auto res = antiMongeDDSMAWK<long long>(N, func, 2);
    cout << res.back().first << endl;
}
void ABC_480_G_CHT() {
    int N;
    cin >> N;
    long long INF = numeric_limits<long long>::max() / 2;
    long long M = numeric_limits<long long>::max() / (N * 3);
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});
    auto func = [&](int l, int r) -> long long { 
        if (r < l + 2) return -M * ((l + 1 - r) * 2 + 1);
        else return stma.get(l, r) - stmi.get(l, r);
    };
    //assert(is_anti_monge(N + 1, N + 1, func));
    
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 2;  // 辺 (i, i + 1) は禁止
    auto res = PSCShortestPathByCHT<long long>(N, func, start);
    cout << res.back().first << endl;
}


int main() {
    //COLOCON_2018_final_C();
    //EDPC_Z("simple_larsch");
    //EDPC_Z("larsch");
    //EDPC_Z("cht");
    //Codeforces_189_C("simple_larsch");
    //Codeforces_189_C("larsch");
    //Codeforces_189_C("cht");
    //yukicoder_705("simple_larsch");
    //yukicoder_705("larsch");
    //yukicoder_705("cht");
    //ABC_218_H();
    //yukicoder_952();
    //Codeforces438_F_enum();
    //Codeforces438_F_alien();
    //ABC_480_G_SMAWK();
    ABC_480_G_CHT();
}