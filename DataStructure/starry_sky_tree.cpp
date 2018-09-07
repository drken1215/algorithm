//
// starry sky tree
//
// verified:
//   AOJ Course Range Query - RMQ and RAQ
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_H&lang=jp
//


/*
    // Construction
    SegTree(n, unity_monoid, unity_lazy)
      ex: starry sky tree (区間加算、区間min取得)


    // Initialization
    init(n): サイズ n に初期化
    set(a, v): a 番目の値を v にセットする
    build(): set した値を元にセグメントツリー全体を構築する、O(n)
    

    // Queries
    update(a, b, v): 区間 [a, b) に v を加算, O(log n)
    get(a, b): 区間 [a, b) の最小値を返す, O(log n)
*/


#include <iostream>
#include <vector>
#include <functional>
using namespace std;


template<class Monoid, class Lazy> struct StarrySky {
    const Monoid INF;
    const Lazy ZERO;
    int SIZE;
    vector<pair<Monoid, int> > dat;
    vector<Lazy> lazy;
    
    StarrySky(int n, const Monoid &inf, const Lazy &zero) : INF(inf), ZERO(zero) {
        init(n);
    }
    void init(int n) {
        SIZE = 1;
        while (SIZE < n) SIZE <<= 1;
        dat.assign(SIZE * 2, {INF, -1});
        lazy.assign(SIZE * 2, ZERO);
    }
    
    /* set, a is 0-indexed */
    void set(int a, const Monoid &v) { dat[a + SIZE] = {v, a}; }
    void build() {
        for (int k = SIZE - 1; k > 0; --k)
            dat[k] = min(dat[k*2], dat[k*2+1]);
    }
    
    /* update [a, b) */
    inline void evaluate(int k) {
        if (lazy[k] == ZERO) return;
        if (k < SIZE) lazy[k*2] += lazy[k], lazy[k*2+1] += lazy[k];
        dat[k].first += lazy[k];
        lazy[k] = ZERO;
    }
    inline void update(int a, int b, const Monoid &v, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b) lazy[k] += v, evaluate(k);
        else if (a < r && l < b) {
            update(a, b, v, k*2, l, (l+r)>>1), update(a, b, v, k*2+1, (l+r)>>1, r);
            dat[k] = min(dat[k*2], dat[k*2+1]);
        }
    }
    inline void update(int a, int b, const Monoid &v) { update(a, b, v, 1, 0, SIZE); }
    
    /* get [a, b) */
    inline pair<Monoid,int> get(int a, int b, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)
            return dat[k];
        else if (a < r && l < b)
            return min(get(a, b, k*2, l, (l+r)>>1), get(a, b, k*2+1, (l+r)>>1, r));
        else
            return {INF, -1};
    }
    inline pair<Monoid,int> get(int a, int b) { return get(a, b, 1, 0, SIZE); }
    inline Monoid operator [] (int a) { return get(a, a+1).first; }
    
    /* debug */
    void print() {
        for (int i = 0; i < SIZE; ++i) { cout << (*this)[i]; if (i != SIZE) cout << ","; }
        cout << endl;
    }
};



int main() {
    int N, Q; scanf("%d %d", &N, &Q);
    StarrySky<long long, long long> sk(N, (1LL<<60), 0);

    // 0 で初期化
    for (int i = 0; i < N; ++i) sk.set(i, 0);
    sk.build();

    // 各クエリ
    for (int q = 0; q < Q; ++q) {
        int type, a, b; scanf("%d %d %d", &type, &a, &b);
        if (type == 0) {
            long long v; scanf("%lld", &v);
            sk.update(a, b+1, v);
        }
        else {
            cout << sk.get(a, b+1).first << endl;
        }
    }
}
