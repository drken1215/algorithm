//
// segment tree
//
// verified:
//   ARC 008 D - タコヤキオイシクナール
//     https://beta.atcoder.jp/contests/arc008/tasks/arc008_4
//
//   ABC 281 E - Least Elements 
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
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


// Segment Tree
template<class Monoid> struct SegTree {
    using Func = function<Monoid(Monoid, Monoid)>;
    int N;
    Func F;
    Monoid IDENTITY;
    int SIZE_R;
    vector<Monoid> dat;

    /* initialization */
    SegTree() {}
    SegTree(int n, const Func f, const Monoid &identity)
    : N(n), F(f), IDENTITY(identity) {
        SIZE_R = 1;
        while (SIZE_R < n) SIZE_R *= 2;
        dat.assign(SIZE_R * 2, IDENTITY);
    }
    void init(int n, const Func f, const Monoid &identity) {  
        N = n;
        F = f;
        IDENTITY = identity;
        SIZE_R = 1;
        while (SIZE_R < n) SIZE_R *= 2;
        dat.assign(SIZE_R * 2, IDENTITY);
    }
    
    /* set, a is 0-indexed */
    /* build(): O(N) */
    void set(int a, const Monoid &v) { dat[a + SIZE_R] = v; }
    void build() {
        for (int k = SIZE_R - 1; k > 0; --k)
            dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    
    /* update a, a is 0-indexed, O(log N) */
    void update(int a, const Monoid &v) {
        int k = a + SIZE_R;
        dat[k] = v;
        while (k >>= 1) dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    
    /* get [a, b), a and b are 0-indexed, O(log N) */
    Monoid get(int a, int b) {
        Monoid vleft = IDENTITY, vright = IDENTITY;
        for (int left = a + SIZE_R, right = b + SIZE_R; left < right; 
        left >>= 1, right >>= 1) {
            if (left & 1) vleft = F(vleft, dat[left++]);
            if (right & 1) vright = F(dat[--right], vright);
        }
        return F(vleft, vright);
    }
    Monoid all_get() { return dat[1]; }
    Monoid operator [] (int a) { return dat[a + SIZE_R]; }
    
    /* get max r that f(get(l, r)) = True (0-indexed), O(log N) */
    /* f(IDENTITY) need to be True */
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += SIZE_R;
        Monoid sum = IDENTITY;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(F(sum, dat[l]))) {
                while (l < SIZE_R) {
                    l = l * 2;
                    if (f(F(sum, dat[l]))) {
                        sum = F(sum, dat[l]);
                        ++l;
                    }
                }
                return l - SIZE_R;
            }
            sum = F(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return N;
    }

    /* get min l that f(get(l, r)) = True (0-indexed), O(log N) */
    /* f(IDENTITY) need to be True */
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = N;
        r += SIZE_R;
        Monoid sum = IDENTITY;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(F(dat[r], sum))) {
                while (r < SIZE_R) {
                    r = r * 2 + 1;
                    if (f(F(dat[r], sum))) {
                        sum = F(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - SIZE_R;
            }
            sum = F(dat[r], sum);
        } while ((r & -r) != r);
        return 0;
    }
    
    /* debug */
    void print() {
        for (int i = 0; i < N; ++i) {
            cout << (*this)[i];
            if (i != N-1) cout << ",";
        }
        cout << endl;
    }
};

////////////////////////////////////////////

// ARC 008 D - タコヤキオイシクナール
void solveARC008D() {
    // 入力
    long long N; 
    int M;
    cin >> N >> M;
    vector<long long> p(M), pls;
    vector<double> a(M), b(M);
    for (int i = 0; i < M; ++i) {
        cin >> p[i] >> a[i] >> b[i]; --p[i];
        pls.push_back(p[i]);
    }
    sort(pls.begin(), pls.end());
    pls.erase(unique(pls.begin(), pls.end()), pls.end());
    
    // セグメント木
    int NN = (int)pls.size();
    auto func = [&](pair<double,double> a, pair<double,double> b) {
        return make_pair(a.first * b.first, a.second * b.first + b.second);
    };
    pair<double, double> identity = make_pair(1, 0);
    SegTree<pair<double, double>> seg(NN, func, identity);

    // 処理
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

void solveABC281E() {
    // 入力
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 座標圧縮を準備
    vector<pair<long long,int>> comp(N);
    for (int i = 0; i < N; ++i) comp[i] = {A[i], i};
    sort(comp.begin(), comp.end());
    vector<int> order(N);
    for (int i = 0; i < N; ++i) order[comp[i].second] = i;

    // セグメント木
    using Monoid = pair<int, long long>;
    auto func = [&](Monoid a, Monoid b) -> Monoid {
        return Monoid(a.first + b.first, a.second + b.second);
    };
    Monoid zero = Monoid(0, 0);
    SegTree<Monoid> seg(N, func, zero);

    // push と pop
    auto push = [&](long long x, int id) -> void {
        seg.update(id, Monoid(1, x));
    };
    auto pop = [&](long long x, int id) -> void {
        seg.update(id, zero);
    };
    auto get = [&]() -> long long {
        int r = seg.max_right([&](Monoid r) { 
            return r.first <= K;
        });
        return seg.get(0, r).second;
    };

    for (int i = 0; i < M; ++i) push(A[i], order[i]);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << get() << " ";
        if (i+M < N) {
            push(A[i+M], order[i+M]); 
            pop(A[i], order[i]);
        }
    }
    cout << endl;
}

int main() {
    //solveARC008D();
    solveABC281E();
}
