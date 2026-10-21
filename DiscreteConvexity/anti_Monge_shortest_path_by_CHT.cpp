//
// 一般化 Convex Hull Trick
//    Monge 性を満たす関数群 f_1, f_2, ..., f_N
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   yukicoder No.705 ゴミ拾い Hard
//     https://yukicoder.me/problems/no/705 
//
// Reference:
//     https://codeforces.com/blog/entry/86731
//     https://yukicoder.me/problems/no/705/editorial
// 


#include <bits/stdc++.h>
using namespace std;


// Convex Hull Trick
/*
    Func：Monge 性を仮定
    - MIN: クエリ x の取りうる最小値
    - MAX: クエリ x の取りうる最大値 (f_i(MAX) がオーバーフローしないように注意）
    - INF: 最大値
 
    - insert (Func f_i): add f_i, O(log N)
    - query (x): min_i{ f_i(x) }, O(log N)
*/
template<class T> struct CHT {
    using Func = function<T(T)>;
    struct Node {
        Func func;
        Node *left, *right;
        Node(const Func& f) : left(nullptr), right(nullptr) {
            func = f;
        }
    };
    
    const T MIN, MAX, INF;
    Node* root;
    
    CHT(T MIN, T MAX, T INF) : MIN(MIN), MAX(MAX), INF(INF), root(nullptr) { }
    Node* insert(Node* p, T low, T high, Func& f) {
        if (!p) return new Node(f);
        if (p->func(low) <= f(low) && p->func(high) <= f(high)) return p;
        if (p->func(low) >= f(low) && p->func(high) >= f(high)) {
            p->func = f;
            return p;
        }
        T mid = (low + high) / 2;
        if (p->func(mid) > f(mid)) swap(p->func, f);
        if (p->func(low) >= f(low)) p->left = insert(p->left, low, mid, f);
        else p->right = insert(p->right, mid, high, f);
        return p;      
    }
    void insert(Func f) {
        root = insert(root, MIN, MAX, f);
    }
    T query(Node* p, T low, T high, T x) {
        if (!p) return INF;
        if (low == high) return p->func(x);
        T mid = (low + high) / 2;
        if (x <= mid) return min(p->func(x), query(p->left, low, mid, x));
        else return min(p->func(x), query(p->right, mid, high, x));
    }
    T query(T x) {
        return query(root, MIN, MAX, x);
    }
};



//------------------------------//
// Examples
//------------------------------//

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    CHT<long long> cht(0, 210000, 1LL<<60);
    for (long long i = 0; i < N; ++i) {
        auto func = [i, &a](long long x) -> long long { return (x - i) * (x - i) + a[i]; };
        cht.insert(func);
    }
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i);
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
    vector<long long> dp(N, INF);
    /*
    　　dp[i] = min_j(dp[j] + (H[j] - H[i])² + C)
    */
    dp[0] = 0;
    CHT<long long> cht(0, MAX, INF);
    auto func = [&H, C](long long x) -> long long {
        return (x - H[0]) * (x - H[0]) + C;
    };
    cht.insert(func);
    for (int i = 1; i < N; i++) {
        long long val = cht.query(H[i]);
        dp[i] = min(dp[i], val);
        auto func = [&dp, &H, i, C](long long x) -> long long {
            return (x - H[i]) * (x - H[i]) + dp[i] + C;
        };
        cht.insert(func);
    }
    long long res = dp[N-1];
    cout << res << endl;
}

// yukicoder No.705 ゴミ拾い Hard
void yukicoder_705() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];

    const long long INF = 1LL<<60;
    const long long MAX = 110000;
    vector<long long> dp(N+1, INF);
    CHT<long long> cht(0, MAX, INF);
    dp[0] = 0;
    for (int i = 1; i <= N; i++) {
        auto func = [&dp, &X, &Y, i](long long x) -> long long {
            long long dx = abs(x - X[i-1]), dy = abs(Y[i-1]);
            return dp[i-1] + dx*dx*dx + dy*dy*dy;
        };
        cht.insert(func);
        dp[i] = cht.query(A[i-1]);
    }
    cout << dp[N] << endl;
}


int main() {
    //COLOCON_2018_final_C();
    //EDPC_Z();
    yukicoder_705();
}