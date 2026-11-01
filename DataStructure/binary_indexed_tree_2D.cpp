//
// Binary Indexed Tree 2D
//
// verified:
//   AOJ 2842 Taiyaki-Master and Eater
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2842
//


#include <iostream>
#include <vector>
#include <queue>
using namespace std;


// BIT 2D
template <class Abel> struct BIT2D {
    const Abel UNITY_SUM = 0;
    vector<vector<Abel> > dat;

    // [0, n) x [0, m)
    BIT2D(int n, int m, Abel unity = 0) : UNITY_SUM(unity),
                                          dat(n, vector<Abel>(m, UNITY_SUM)) { }
    void init(int n, int m) {
        dat.assign(n, vector<Abel>(m, UNITY_SUM));
    }
    
    // add x on the point (a, b)
    inline void add(int a, int b, Abel x) {
        for (int i = a; i < (int)dat.size(); i |= i + 1)
            for (int j = b; j < (int)dat[0].size(); j |= j + 1)
                dat[i][j] = dat[i][j] + x;
    }

    // [0, p) x [0, q), 0-indexed
    inline Abel sum(int p, int q) {
        Abel res = UNITY_SUM;
        for (int i = p - 1; i >= 0; i = (i & (i + 1)) - 1)
            for (int j = q - 1; j >= 0; j = (j & (j + 1)) - 1)
                res = res + dat[i][j];
        return res;
    }
    
    // x1 <= x < x2, y1 <= y < y2, 0-indexed
    inline Abel sum(int x1, int x2, int y1, int y2) {
        return sum(x2, y2) - sum(x1, y2) - sum(x2, y1) + sum(x1, y1);
    }
    
    // debug
    void print() {
        for (int i = 1; i < (int)dat.size(); ++i) {
            for (int j = 1; j < (int)dat[0].size(); ++j)
                cout << sum(i, j, i + 1, j + 1) << ",";
            cout << endl;
        }
    }
};



//------------------------------//
// Examples
//------------------------------//

typedef pair<int,int> pint; // zahyou
typedef pair<int,pint> ppint; // (time, zahyou)

int main() {
    int H, W, T, Q;
    cin >> H >> W >> T >> Q;
    BIT2D<int> fin(H, W), mada(H, W);
    queue<ppint> que;
    for (int query = 0; query < Q; ++query) {
        int t, c, h, w;
        cin >> t >> c >> h >> w;
        --h, --w;
        while (!que.empty()) {
            ppint pall = que.front();
            if (pall.first > t) break;
            que.pop();
            int x = pall.second.first, y = pall.second.second;
            mada.add(x, y, -1);
            fin.add(x, y, 1);
        }
        if (c == 0) {
            mada.add(h, w, 1);
            que.push(ppint(t+T, pint(h, w)));
        }
        else if (c == 1) {
            if (fin.sum(h, h+1, w, w+1) == 1) fin.add(h, w, -1);
        }
        else {
            int h2, w2;
            cin >> h2 >> w2;
            --h2, --w2;
            cout << fin.sum(h, h2+1, w, w2+1) << " "
                 << mada.sum(h, h2+1, w, w2+1) << endl;
        }
    }
}
