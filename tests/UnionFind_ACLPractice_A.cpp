#define PROBLEM "https://atcoder.jp/contests/practice2/tasks/practice2_a"
#define IGNORE

#include <bits/stdc++.h>
using namespace std;

#include "DataStructure/UnionFind.hpp"
    
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
        //cout << uf << endl;
    }
}