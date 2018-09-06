//
// Binary Indexed Tree 2D
//
// verified:
//   AOJ 2842 たい焼きマスターと食べ盛り
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2842
//


#include <iostream>
#include <vector>
#include <queue>
using namespace std;


template <class Abel> struct BIT2D {
    const Abel UNITY_SUM = 0;                       // to be set
    vector<vector<Abel> > dat;
    
    BIT2D(int n = 1, int m = 1) : dat(n + 1, vector<Abel>(m + 1, UNITY_SUM)) { }
    void init(int n, int m) { dat.assign(n + 1, vector<Abel>(m + 1, UNITY_SUM)); }
    
    /* add x on the point (a, b) */
    inline void add(int a, int b, Abel x) {
        for (int i = a; i < (int)dat.size(); i += i & -i)
            for (int j = b; j < (int)dat[0].size(); j += j & -j)
                dat[i][j] = dat[i][j] + x;
    }
    
    inline Abel sum(int p, int q) {
        Abel res = UNITY_SUM;
        for (int i = p; i > 0; i -= i & -i)
            for (int j = q; j > 0; j -= j & -j)
                res = res + dat[i][j];
        return res;
    }
    
    /* a1 <= x < b1, a2 <= y < b2, 1-indexd */
    inline Abel sum(int a1, int a2, int b1, int b2) {
        return sum(b1-1, b2-1) - sum(a1-1, b2-1) - sum(b1-1, a2-1) + sum(a1-1, a2-1);
    }
    
    /* debug */
    void print() {
        for (int i = 1; i < (int)dat.size(); ++i) {
            for (int j = 1; j < (int)dat[0].size(); ++j)
                cout << sum(i, j, i+1, j+1) << ",";
            cout << endl;
        }
    }
};



typedef pair<int,int> pint; // zahyou
typedef pair<int,pint> ppint; // (time, zahyou)

int main() {
    int H, W, T, Q; cin >> H >> W >> T >> Q;
    BIT2D<int> fin(H, W), mada(H, W);
    queue<ppint> que;
    for (int query = 0; query < Q; ++query) {
        int t, c, h, w; cin >> t >> c >> h >> w;
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
            if (fin.sum(h, w, h+1, w+1) == 1) fin.add(h, w, -1);
        }
        else {
            int h2, w2; cin >> h2 >> w2;
            cout << fin.sum(h, w, h2+1, w2+1) << " " << mada.sum(h, w, h2+1, w2+1) << endl;
        }
    }
}
