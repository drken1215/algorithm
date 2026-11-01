//
// Disjoint Sparse Table
//
// cf.
//   noshi: Disjoint Sparse Table と セグ木に関するポエム
//     http://noshi91.hatenablog.com/entry/2018/05/08/183946
//
// verified
//   Codechef Product on the segment by modulo
//     https://www.codechef.com/problems/SEGPROD
//


/*
  Sparse Table は本来は「交叉半束」上で動作
  Disjoint Sparse Table は「半群」上でも動作できるようにしたバージョン
    静的な半群列の区間和を高速に計算する
 
  半群: 二項演算が定義されて、結合法則を満たす、単位元も逆元も不要 (min, max, 和, 積, gcd)
   (モノイド: 半群 + 単位元)
   (群: 半群 + 単位元 + 逆元)
  交叉半束: さらに「べき等則 (a + a = a)」「可換則 (a + b = b + a)」を満たす必要がある

  init(vec): 配列を vec で初期化構築, O(n logn)
  get(a, b): [a, b) の区間和を取得
*/


#include <iostream>
#include <vector>
#include <functional>
using namespace std;


template<class SemiGroup> struct DisjointSparseTable {
    using Func = function< SemiGroup(SemiGroup, SemiGroup) >;
    const Func F;
    vector<vector<SemiGroup> > dat;
    vector<int> height;
    
    DisjointSparseTable(const Func &f) : F(f) { }
    DisjointSparseTable(const Func &f, const vector<SemiGroup> &vec) : F(f) { init(vec); }
    void init(const vector<SemiGroup> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<SemiGroup>(n));
        height.assign((1<<h), 0);
        for (int i = 2; i < (1<<h); i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            int s = (1<<i);
            for (int j = 0; j < n; j += (s << 1)) {
                int t = min(j+s, n);
                dat[i][t-1] = vec[t-1];
                for (int k = t-2; k >= j; --k) dat[i][k] = F(vec[k], dat[i][k+1]);
                if (n <= t) break;
                dat[i][t] = vec[t];
                for (int k = t+1; k < min(t+s, n); ++k) dat[i][k] = F(dat[i][k-1], vec[k]);
            }
        }
    }
    SemiGroup get(int a, int b) {
        if (a >= --b) return dat[0][a];
        return F(dat[height[a^b]][a], dat[height[a^b]][b]);
    }
};



//------------------------------//
// Examples
//------------------------------//

int main() {
    int T; scanf("%d", &T);
    int N, P, Q;
    DisjointSparseTable<long long> dst([&](long long a, long long b){return a*b%P;});
    for (int CASE = 0; CASE < T; ++CASE) {
        scanf("%d %d %d", &N, &P, &Q);
        vector<long long> v(N);
        for (int i = 0; i < N; ++i) scanf("%lld", &v[i]);
        dst.init(v);
        vector<int> b(Q/64+2);
        for (int i = 0; i < (int)b.size(); ++i) scanf("%d", &b[i]);
        int x = 0, L = 0, R = 0;
        for (int query = 0; query < Q; ++query) {
            if (query % 64 == 0) L = (b[query/64] + x) % N, R = (b[query/64+1] + x) % N;
            else L = (L + x) % N, R = (R + x) % N;
            if (L > R) swap(L, R);
            x = (dst.get(L, R+1) + 1) % P;
        }
        printf("%d\n", x);
    }
}
