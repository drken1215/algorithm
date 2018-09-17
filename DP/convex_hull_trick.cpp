//
// Convex Hull Trick
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//


#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;


// Convex Hull Trick
/*
    直線の傾きに単調性を仮定しない
    - MIN: クエリ x の取りうる最小値
    - MAX: クエリ x の取りうる最大値
    - INF: 最大値
 
    - insert (a, b): add y = ax + b
    - query (x): min_i{ l_i->get(x) }
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
        Line l;
        Node *lhs,*rhs;
        Node(Line l) : l(l), lhs(nullptr), rhs(nullptr){}
    };
    
    const T MIN, MAX, INF;
    Node* root;
    
    CHT(T MIN, T MAX, T INF) : MIN(MIN), MAX(MAX), INF(INF), root(nullptr) { }
    Node* insert(Node* p, T left, T right, Line& l){
        if(!p) return new Node(l);
        if (p->l.get(left) <= l.get(left) && p->l.get(right) <= l.get(right)) return p;
        if (p->l.get(left) >= l.get(left) && p->l.get(right) >= l.get(right)) {
            p->l = l;
            return p;
        }
        T mid = (left + right) / 2;
        if (p->l.get(mid) > l.get(mid)) swap(p->l , l);
        if (p->l.get(left) >= l.get(left)) p->lhs = insert(p->lhs, left, mid, l);
        else p->rhs = insert(p->rhs, mid+1, right, l);
        return p;
        
    }
    void insert(T a, T b){
        Line l(a, b);
        root = insert(root, MIN, MAX, l);
    }
    T query(Node* p, T left, T right, T x){
        if(!p) return INF;
        if(left == right) return p->l.get(x);
        T mid = (left + right) / 2;
        if (x <= mid) return min(p->l.get(x), query(p->lhs, left, mid, x));
        else return min(p->l.get(x), query(p->rhs, mid+1, right, x));
    }
    T query(T x){
        return query(root, MIN, MAX, x);
    }
};


int main() {
    long long N; cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    CHT<long long> cht(0, 210000, 1LL<<60);
    for (long long i = 0; i < N; ++i) cht.insert(-2LL*i, a[i] + i*i);
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i) + i*i;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}
