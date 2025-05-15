//
// Convex Hull Trick (追加クエリの直線の傾きが単調減少の場合)
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   yukicoder No.952 危険な火薬庫
//     https://yukicoder.me/problems/no/952
//


#include <bits/stdc++.h>
using namespace std;


// Convex Hull Trick
/*
    追加クエリの直線の傾きが単調減少の場合
 
    - insert(a, b): add y = ax + b, O(1)
    - query(x): min_i{a[i]*x + b}, O(log N)
    - query_monotone(x): min_i{a[i]*x + b}, O(1), x も単調増加であることを仮定
*/
template<class T> struct CHT {
    using Line = pair<T, T>;
    deque<Line> lines;
    T INF;

    CHT(T inf = numeric_limits<T>::max() / 2) : INF(inf) {}
    bool check(const Line &a, const Line &b, const Line &c) {
        return (b.first - a.first) * (c.second - b.second)
        >= (b.second - a.second) * (c.first - b.first);
    }
    
    T get(int k, T x) {
        if (lines.empty() || k >= lines.size()) return INF;
        assert(k >= 0 && k < lines.size());
        return lines[k].first * x + lines[k].second;
    }
    
    void insert(T a, T b) {
        Line l(a, b);
        while (lines.size() >= 2
               && check(lines[(int)lines.size()-2], lines[(int)lines.size()-1], l)) 
            lines.pop_back();
        lines.push_back(l);
    }
    
    T query(T x) {
        int low = -1, high = (int)lines.size();
        while (high - low > 1) {
            int mid = (low + high) / 2;
            if (get(mid, x) >= get(mid + 1, x)) low = mid;
            else high = mid;
        }
        return get(high, x);
    }
    
    // クエリの単調性も成り立つ場合 (x が単調増加)
    T query_monotone(T x) {
        while (lines.size() >= 2 && get(0, x) >= get(1, x)) lines.pop_front();
        if (lines.empty()) return INF;
        return lines[0].first * x + lines[0].second;
    }
};



//------------------------------//
// Examples
//------------------------------//

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N, 1LL<<60);
    for (int i = 0; i < N; ++i) cin >> a[i];
    CHT<long long> cht;
    for (long long i = 0; i < N; ++i) cht.insert(-2LL * i, a[i] + i*i);
    for (long long i = 0; i < N; ++i) res[i] = cht.query_monotone(i) + i*i;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}

// AtCoder EDPC Z - Frog 3
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (int i = 0; i < N; i++) cin >> H[i];

    const long long INF = 1LL<<60;
    const long long MAX = 1LL<<40;
    vector<long long> dp2(N, INF);
    /*
    　　dp[i] = min_j(dp[j] + (H[j] - H[i])² + C)
    　　dp2[i] = dp[i] + H[i]² とすると
    　　dp2[i] = min_j(-2 H[j] × H[i] + dp2[j]) + 2 H[i]² + C
    */
    dp2[0] = H[0] * H[0];
    CHT<long long> cht(INF);
    cht.insert(-H[0] * 2, dp2[0]);
    for (int i = 1; i < N; i++) {
        long long val = cht.query_monotone(H[i]);
        dp2[i] = min(dp2[i], val + H[i] * H[i] * 2 + C);
        cht.insert(-H[i] * 2, dp2[i]);
    }
    long long res = dp2[N-1] - H[N-1] * H[N-1];
    cout << res << endl;
}

// yukicoder No.952 危険な火薬庫
using i128 = __int128;
constexpr i128 to_integer(const string &s) {
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
void yukicoder_952() {
    int N;
    cin >> N;
    N++;
    vector<i128> A(N, 0), S(N+1, 0);
    for (int i = 0; i < N; i++) {
        if (i) cin >> A[i];
        S[i+1] = S[i] + A[i];
    }

    vector chts(N+1, CHT<i128>());
    
    i128 INF = i128(1)<<100;
    vector dp(N+1, vector(N+1, INF));

    auto push = [&](int i, int j) -> void {
        if (i-1 >= 0) chts[i-j].insert(-S[i]*2, S[i]*S[i] + dp[i-1][j]);
    };

    dp[0][0] = 0;
    push(0, 0);

    for (int i = 1; i <= N; i++) {
        for (int j = 0; j < i; j++) {
            dp[i][j] = min(dp[i][j], dp[i-1][j]);
            dp[i][j] = min(dp[i][j], chts[i-j].query_monotone(S[i]) + S[i] * S[i]);
            push(i, j);
        }
    }

    for (int k = 1; k < N; k++) {
        cout << dp[N][k] << endl;
    }
}


int main() {
    //COLOCON_2018_final_C();
    //EDPC_Z();
    yukicoder_952();
}