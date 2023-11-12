//
// next combination (nCk 通りの全探索)
//
// reference:
//   ビット演算 (bit 演算) の使い方を総特集！ 〜 マスクビットから bit DP まで 〜
//     https://qiita.com/drken/items/7c6ff2aa4d8fce1c9361
//
// verifed:
//   AtCoder ABC 328 E - Modulo MST
//     https://atcoder.jp/contests/abc328/tasks/abc328_e
//


#include <bits/stdc++.h>
using namespace std;


// next combination
template<class T> bool next_combination(T &bit, int N) {
    T x = bit & -bit, y = bit + x;
    bit = (((bit & ~y) / x) >> 1) | y;
    return (bit < (1LL << N));
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// small test
void small_test() {
    int N = 5;  // {0, 1, 2, 3, 4} の部分集合を考える
    int K = 3;
    
    int bit = (1 << K) - 1;  // bit = {0, 1, 2}
    do {
        cout << bitset<5>(bit) << endl;
    } while (next_combination(bit, N));
}


// ABC 328 E
void ABC_328_E() {
    // edge
    struct Edge {
        int from, to;
        long long weight;
        Edge(int f = 0, int t = 0, long long w = 0) : from(f), to(t), weight(w) {}
    };
    struct UnionFind {
        // core member
        vector<int> par;

        // constructor
        UnionFind() { }
        UnionFind(int n) : par(n, -1) { }
        void init(int n) { par.assign(n, -1); }
        
        // core methods
        int root(int x) {
            if (par[x] < 0) return x;
            else return par[x] = root(par[x]);
        }
        
        bool same(int x, int y) {
            return root(x) == root(y);
        }
        
        bool merge(int x, int y) {
            x = root(x), y = root(y);
            if (x == y) return false;
            if (par[x] > par[y]) swap(x, y); // merge technique
            par[x] += par[y];
            par[y] = x;
            return true;
        }
        
        int size(int x) {
            return -par[root(x)];
        }
    };
    
    long long N, M, K;
    cin >> N >> M >> K;
    vector<Edge> edges(M);
    for (int i = 0; i < M; ++i) {
        long long u, v, w;
        cin >> u >> v >> w;
        --u, --v;
        edges[i] = Edge(u, v, w);
    }
    
    long long res = 1LL<<62;
    long long bit = (1LL << (N-1)) - 1;
    do {
        long long tmp = 0;
        bool ok = true;
        UnionFind uf(N);
        for (int i = 0; i < M; ++i) {
            if (bit >> i & 1) {
                auto [f, t, w] = edges[i];
                if (uf.same(f, t)) ok = false;
                uf.merge(f, t);
                tmp += w;
            }
        }
        if (ok) res = min(res, tmp % K);
    } while (next_combination(bit, M));
    
    cout << res << endl;
}


int main() {
    //small_test();
    ABC_328_E();
}

