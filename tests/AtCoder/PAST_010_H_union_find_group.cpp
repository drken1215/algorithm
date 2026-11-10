#define PROBLEM "https://atcoder.jp/contests/past202203-open/tasks/past202203_h"

#include "union_find.hpp"

int main() {
    int N, Q, x, y;
    cin >> N >> Q;
    UnionFind uf(N);
    for (int i = 0; i < Q; i++) {
        int typ;
        cin >> typ;
        if (typ == 1) {
            cin >> x >> y, x--, y--;
            uf.merge(x, y);
        } else {
            cin >> x, x--;
            auto res = uf.group(x);
            sort(res.begin(), res.end());
            for (auto v : res) cout << v+1 << " ";
            cout << endl;
        }
    }
}
