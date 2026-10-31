//
// Union-Find (union by size)
//
// verified:
//   Yosupo Library Checker - Unionfind
//     https://judge.yosupo.jp/problem/unionfind
//
//   AOJ Course DSL_1_A Set - Disjoint Set: Union Find Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_1_A&lang=jp
//


/*
    same(x, y): x と y が同じグループにいるか, 計算量はならし O(α(n))
    merge(x, y): x と y を同じグループにする, 計算量はならし O(α(n))
    size(x): x を含むグループの所属メンバー数, 計算量はならし O(α(n))
    groups(): グループ分けの構造を返す, 計算量は O(n)
*/


#include <bits/stdc++.h>
using namespace std;

#include "../DataStructureBasic/union_find_size.hpp"

//------------------------------//
// Examples
//------------------------------//
    
int main() {
    int N, Q;
    cin >> N >> Q;
    UnionFind uf(N);
    for (int q = 0; q < Q; ++q) {
        int type, x, y;
        cin >> type >> x >> y;
        if (type == 0) uf.merge(x, y);
        else {
            if (uf.same(x, y)) cout << 1 << endl;
            else cout << 0 << endl;
        }
        cout << uf << endl;
    }
}