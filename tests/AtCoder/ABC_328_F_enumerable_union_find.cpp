#define PROBLEM "https://atcoder.jp/contests/abc328/tasks/abc328_f"

#include "union_find.hpp"

int main() {
    int N, Q;
    cin >> N >> Q;
    UnionFind uf(N);
    vector<long long> x(N, 0);  // x value
    for (int i = 0; i < Q; ++i) {
        int a, b, d;a
        cin >> a >> b >> d;
        --a, --b;
        
        bool good = true;
        if (uf.same(a, b)) {
            // x[a] - x[b] = d でないとき、ダメ
            if (x[a] - x[b] != d) good = false;
        } else {
            // マージテクにより、size(a) < size(b) となるようにする
            if (uf.size(a) > uf.size(b)) {
                swap(a, b);
                d = -d;
            }
            
            // a を含むグループのメンバーたちに足す値を求める
            long long add = (x[b] + d) - x[a];
            
            // x の値を調整して、マージする
            vector<int> ids = uf.group(a);
            for (int id : ids) x[id] += add;
            uf.merge(a, b);
        }

        // good ならば出力する
        if (good) cout << i+1 << " ";
    }
    cout << endl;
}
