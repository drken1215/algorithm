//
// Convex Hull Trick
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//


#include <bits/stdc++.h>
using namespace std;


// Convex Hull Trick
/*
    直線の傾きに単調性を仮定しない
    - MIN: クエリ x の取りうる最小値
    - MAX: クエリ x の取りうる最大値 (a・MAX + b がオーバーフローしないように注意）
    - INF: 最大値
 
    - insert (a, b): add y = ax + b, O(log N)
    - query (x): min_i{ l_i->get(x) }, O(log N)
*/
template<class T> struct CHT {
    struct Line {
        T a, b;
        Line(T a = 0,T b = 0) : a(a), b(b) { }
        T get(T x) {
            return a*x + b;
        }
    };
    
    struct Node {
        Line line;
        Node *left, *right;
        Node(Line line) : line(line), left(nullptr), right(nullptr){}
    };
    
    const T MIN, MAX, INF;
    Node* root;
    
    CHT(T MIN, T MAX, T INF) : MIN(MIN), MAX(MAX), INF(INF), root(nullptr) { }
    Node* insert(Node* p, T low, T high, Line& l){
        if(!p) return new Node(l);
        if (p->line.get(low) <= l.get(low) && p->line.get(high) <= l.get(high)) return p;
        if (p->line.get(low) >= l.get(low) && p->line.get(high) >= l.get(high)) {
            p->line = l;
            return p;
        }
        T mid = (low + high) / 2;
        if (p->line.get(mid) > l.get(mid)) swap(p->line , l);
        if (p->line.get(low) >= l.get(low)) p->left = insert(p->left, low, mid, l);
        else p->right = insert(p->right, mid, high, l);
        return p;      
    }
    void insert(T a, T b){
        Line l(a, b);
        root = insert(root, MIN, MAX, l);
    }
    T query(Node* p, T low, T high, T x){
        if(!p) return INF;
        if(low == high) return p->line.get(x);
        T mid = (low + high) / 2;
        if (x <= mid) return min(p->line.get(x), query(p->left, low, mid, x));
        else return min(p->line.get(x), query(p->right, mid, high, x));
    }
    T query(T x){
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
    for (long long i = 0; i < N; ++i) cht.insert(-2LL*i, a[i] + i*i);
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i) + i*i;
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
    CHT<long long> cht(0, MAX, INF);
    cht.insert(-H[0] * 2, dp2[0]);
    for (int i = 1; i < N; i++) {
        long long val = cht.query(H[i]);
        dp2[i] = min(dp2[i], val + H[i] * H[i] * 2 + C);
        cht.insert(-H[i] * 2, dp2[i]);
    }
    long long res = dp2[N-1] - H[N-1] * H[N-1];
    cout << res << endl;
}


int main() {
    //COLOCON_2018_final_C();
    EDPC_Z();
}