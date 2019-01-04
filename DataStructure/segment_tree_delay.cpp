//
// segment tree (with delay)
//
// verified:
//   AOJ Course Range Query - RMQ and RAQ
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_H&lang=jp
//


/*
    セグメントツリーは「作用つきモノイド」上で定義される
      FM(a, b): 2 つのモノイド間に定義される演算
      FA(a, d): モノイド元 a への作用素 d による作用
      FL(d, e): 作用素 d への作用素 e による作用
      UNITY_MONOID: モノイドの単位元
      UNITY_LAZY: 作用素の単位元


    // Construction
    SegTree(n, fm, fa, fl, unity_monoid, unity_lazy)
      ex: starry sky tree (区間加算、区間min取得)
        auto fm = [](long long a, long long b) { return min(a, b); };
        auto fa = [](long long &a, long long d) { a += d; };
        auto fl = [](long long &d, long long e) { d += e; };
        SegTree<long long, long long> seg(N, fm, fa, fl, (1LL<<60), 0);


    // Initialization
    init(n): サイズ n に初期化
    set(a, v): a 番目の値を v にセットする
    build(): set した値を元にセグメントツリー全体を構築する、O(n)
    

    // Queries
    update(a, b, v): 区間 [a, b) を作用素 v を用いて更新する, O(log n)
    get(a, b): 区間 [a, b) についての演算結果を返す, O(log n)
*/


#include <iostream>
#include <vector>
#include <functional>
using namespace std;


// Segment Tree
template<class Monoid, class Action> struct SegTree {
    using FuncMonoid = function< Monoid(Monoid, Monoid) >;
    using FuncAction = function< void(Monoid&, Action) >;
    using FuncLazy = function< void(Action&, Action) >;
    FuncMonoid FM;
    FuncAction FA;
    FuncLazy FL;
    Monoid UNITY_MONOID;
    Action UNITY_LAZY;
    int SIZE, HEIGHT;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    SegTree() { }
    SegTree(int n, const FuncMonoid fm, const FuncAction fa, const FuncLazy fl,
            const Monoid &unity_monoid, const Action &unity_lazy)
    : FM(fm), FA(fa), FL(fl), UNITY_MONOID(unity_monoid), UNITY_LAZY(unity_lazy) {
        SIZE = 1; HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, UNITY_MONOID);
        lazy.assign(SIZE * 2, UNITY_LAZY);
    }
    void init(int n, const FuncMonoid fm, const FuncAction fa, const FuncLazy fl,
              const Monoid &unity_monoid, const Action &unity_lazy) {
        FM = fm; FA = fa; FL = fl;
        UNITY_MONOID = unity_monoid; UNITY_LAZY = unity_lazy;
        SIZE = 1; HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, UNITY_MONOID);
        lazy.assign(SIZE * 2, UNITY_LAZY);
    }
    
    /* set, a is 0-indexed */
    void set(int a, const Monoid &v) { dat[a + SIZE] = v; }
    void build() {
        for (int k = SIZE - 1; k > 0; --k)
            dat[k] = FM(dat[k*2], dat[k*2+1]);
    }
    
    /* update [a, b) */
    inline void evaluate(int k) {
        if (lazy[k] == UNITY_LAZY) return;
        if (k < SIZE) FL(lazy[k*2], lazy[k]), FL(lazy[k*2+1], lazy[k]);
        FA(dat[k], lazy[k]);
        lazy[k] = UNITY_LAZY;
    }
    inline void update(int a, int b, const Action &v, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)  FL(lazy[k], v), evaluate(k);
        else if (a < r && l < b) {
            update(a, b, v, k*2, l, (l+r)>>1), update(a, b, v, k*2+1, (l+r)>>1, r);
            dat[k] = FM(dat[k*2], dat[k*2+1]);
        }
    }
    inline void update(int a, int b, const Action &v) { update(a, b, v, 1, 0, SIZE); }
    
    /* get [a, b) */
    inline Monoid get(int a, int b, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)
            return dat[k];
        else if (a < r && l < b)
            return FM(get(a, b, k*2, l, (l+r)>>1), get(a, b, k*2+1, (l+r)>>1, r));
        else
            return UNITY_MONOID;
    }
    inline Monoid get(int a, int b) { return get(a, b, 1, 0, SIZE); }
    inline Monoid operator [] (int a) { return get(a, a+1); }
    
    /* debug */
    void print() {
        for (int i = 0; i < SIZE; ++i) { cout << (*this)[i]; if (i != SIZE) cout << ","; }
        cout << endl;
    }
};


int main() {
    int N, Q; scanf("%d %d", &N, &Q);
    auto fm = [](long long a, long long b) { return min(a, b); };
    auto fa = [](long long &a, long long d) { a += d; };
    auto fl = [](long long &d, long long e) { d += e; };
    SegTree<long long, long long> seg(N, fm, fa, fl, (1LL<<60), 0);

    // 0 で初期化
    for (int i = 0; i < N; ++i) seg.set(i, 0);
    seg.build();

    // 各クエリ
    for (int q = 0; q < Q; ++q) {
        int type, a, b; scanf("%d %d %d", &type, &a, &b);
        if (type == 0) {
            long long v; scanf("%lld", &v);
            seg.update(a, b+1, v);
        }
        else {
            cout << seg.get(a, b+1) << endl;
        }
        //seg.print();
    }
}
