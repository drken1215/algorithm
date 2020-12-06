//
// segment tree
//
// verified:
//   ARC 008 D - タコヤキオイシクナール
//     https://beta.atcoder.jp/contests/arc008/tasks/arc008_4
//


/*
    セグメントツリーは二項演算の定義されたモノイド上で定義される
      二項演算関数 f(x, y) を構造体に渡す


    // Construction
    SegTree(n, f, unity): サイズ n に初期化、f は二項演算、unity は単位元 (min なら INF, + なら 0) 
      ex
      ・区間和: SegTree<int> seg(n, [](int a, int b){ return a + b; }, 0);
      ・区間min: SegTree<int> seg(n, [](int a, int b}{ return min(a, b); }, INF);


    // Initialization
    init(n): サイズ n に初期化
    set(a, v): a 番目の値を v にセットする
    build(): set した値を元にセグメントツリー全体を構築する、O(n)
    

    // Queries
    update(a, v): a 番目の値を v に更新する, O(log n)
    get(a, b): 区間 [a, b) についての演算結果を返す, O(log n)
*/


#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <iomanip>
using namespace std;


template<class Monoid> struct SegTree {
    using Func = function<Monoid(Monoid, Monoid)>;
    const Func F;
    const Monoid UNITY;
    int SIZE_R;
    vector<Monoid> dat;

    SegTree() {}
    SegTree(int n, const Func f, const Monoid &unity): F(f), UNITY(unity) { init(n); }
    void init(int n) {
        SIZE_R = 1;
        while (SIZE_R < n) SIZE_R *= 2;
        dat.assign(SIZE_R * 2, UNITY);
    }
    
    /* set, a is 0-indexed */
    void set(int a, const Monoid &v) { dat[a + SIZE_R] = v; }
    void build() {
        for (int k = SIZE_R - 1; k > 0; --k)
            dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    
    /* update a, a is 0-indexed */
    void update(int a, const Monoid &v) {
        int k = a + SIZE_R;
        dat[k] = v;
        while (k >>= 1) dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    
    /* get [a, b), a and b are 0-indexed */
    Monoid get(int a, int b) {
        Monoid vleft = UNITY, vright = UNITY;
        for (int left = a + SIZE_R, right = b + SIZE_R; left < right; left >>= 1, right >>= 1) {
            if (left & 1) vleft = F(vleft, dat[left++]);
            if (right & 1) vright = F(dat[--right], vright);
        }                                                                                                              
        return F(vleft, vright);
    }
    inline Monoid operator [] (int a) { return dat[a + SIZE_R]; }
    
    /* debug */
    void print() {
        for (int i = 0; i < SIZE_R; ++i) {
            cout << (*this)[i];
            if (i != SIZE_R-1) cout << ",";
        }
        cout << endl;
    }
};



int main() {
    long long N; int M; cin >> N >> M;
    vector<long long> p(M), pls;
    vector<double> a(M), b(M);
    for (int i = 0; i < M; ++i) {
        cin >> p[i] >> a[i] >> b[i]; --p[i];
        pls.push_back(p[i]);
    }
    sort(pls.begin(), pls.end());
    pls.erase(unique(pls.begin(), pls.end()), pls.end());
    int NN = (int)pls.size();

    SegTree<pair<double,double> > seg(NN,
                                      [](pair<double,double> a, pair<double,double> b){
                                          return make_pair(a.first * b.first, a.second * b.first + b.second);
                                      },
                                      make_pair(1, 0)
                                      );

    double Min = 1.0, Max = 1.0;
    for (int i = 0; i < M; ++i) {
        int idx = lower_bound(pls.begin(), pls.end(), p[i]) - pls.begin();
        seg.update(idx, make_pair(a[i], b[i]));
        pair<double,double> res = seg.get(0, NN);
        Min = min(Min, res.first + res.second);
        Max = max(Max, res.first + res.second);
    }
    cout << fixed << setprecision(10) << Min << endl << Max << endl;
}
