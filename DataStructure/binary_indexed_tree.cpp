//
// Binary Indexed Tree
//
// verified:
//   AOJ Course DSL_2_B Range Query - Range Sum Query (RSQ)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_B&lang=jp
//

#include <iostream>
#include <vector>
using namespace std;


template <class Abel> struct BIT {
    vector<Abel> dat;
    Abel UNITY_SUM = 0;						// to be set
    
    /* [1, n] */
    BIT(int n) : dat(n + 1, UNITY_SUM) { }
    void init(int n) { dat.sssign(n + 1, UNITY_SUM); }
    
    /* a is 1-indexed */
    inline void add(int a, Abel x) {
        for (int i = a; i < (int)dat.size(); i += i & -i)
            dat[i] = dat[i] + x;
    }
    
    /* [1, a], a is 1-indexed */
    inline Abel sum(int a) {
        Abel res = UNITY_SUM;
        for (int i = a; i > 0; i -= i & -i)
            res = res + dat[i];
        return res;
    }

    /* [a, b), a and b are 1-indexed */
    inline Abel sum(int a, int b) {
        return sum(b - 1) - sum(a - 1);
    }
    
    /* debug */
    void print() {
        for (int i = 1; i < (int)dat.size(); ++i) cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};



int main() {
    int N, Q; cin >> N >> Q;
    BIT<int> bit(N);
    for (int query = 0; query < Q; ++query) {
        int com, x, y; cin >> com >> x >> y;
        if (com == 0) bit.add(x, y);
        else cout << bit.sum(x, y+1) << endl;
    }
}
