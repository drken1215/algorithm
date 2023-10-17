//
// RMQ (by segment tree)
//
// verified:
//   AOJ Course DSL_2_A Range Query - Range Minimum Query (RMQ)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_A&lang=jp
//

/*
    RMQ(n, inf): サイズ n に初期化、無限大の値を inf とする
    init(n): サイズ n に初期化
    set(a, v): a 番目の値を v にセットする
    build(): set した値を元にセグメントツリー全体を構築する、O(n)
    
    update(a, v): a 番目の値を v に更新する, O(log n)
    get(a, b): 区間 [a, b) の (最小値, 最小値をとる index) を返す, O(log n)
*/


#include <iostream>
#include <vector>
using namespace std;


template<class Monoid> struct RMQ {
    Monoid INF;
    int SIZE_R;
    vector<pair<Monoid,int> > dat;
    
    RMQ() {}
    RMQ(int n, const Monoid &inf): INF(inf) { 
        init(n, inf);
    }
    void init(int n, const Monoid &inf) {
        INF = inf;
        SIZE_R = 1;
        while (SIZE_R < n) SIZE_R *= 2;
        dat.assign(SIZE_R * 2, pair<Monoid,int>(INF, -1));
    }
    
    /* set, a is 0-indexed */
    void set(int a, const Monoid &v) { dat[a + SIZE_R] = make_pair(v, a); }
    void build() {
        for (int k = SIZE_R - 1; k > 0; --k) {
            dat[k] = min(dat[k*2], dat[k*2+1]);
        }
    }
    
    /* update, a is 0-indexed */
    void update(int a, const Monoid &v) {
        int k = a + SIZE_R;
        dat[k] = make_pair(v, a);
        while (k >>= 1) dat[k] = min(dat[k*2], dat[k*2+1]);
    }
    
    /* get {min-value, min-index}, a and b are 0-indexed */
    pair<Monoid,int> get(int a, int b) {
        pair<Monoid,int> vleft = make_pair(INF, -1), vright = make_pair(INF, -1);
        for (int left = a + SIZE_R, right = b + SIZE_R; left < right; left >>= 1, right >>= 1) {
            if (left & 1) vleft = min(vleft, dat[left++]);
            if (right & 1) vright = min(dat[--right], vright);
        }
        return min(vleft, vright);
    }
    inline Monoid operator [] (int a) { return dat[a + SIZE_R].first; }
    
    /* debug */
    void print() {
        for (int i = 0; i < SIZE_R; ++i) {
            Monoid val = (*this)[i];
            if (val < INF) cout << val;
            else cout << "INF";
            if (i != SIZE_R-1) cout << ",";
        }
        cout << endl;
    }
};



int main() {
    int N, Q; cin >> N >> Q;
    RMQ<int> rmq(N, 2147483647);
    for (int query = 0; query < Q; ++query) {
        int com, x, y; cin >> com >> x >> y;
        if (com == 0) rmq.update(x, y);
        else cout << rmq.get(x, y+1).first << endl;
    }
}
