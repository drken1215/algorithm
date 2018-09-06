//
// Binary Indexed Tree (「区間加算」「区間和取得」両対応)
//
// verified:
//   AOJ Course DSL_2_B Range Query - Range Sum Query (RSQ)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_B&lang=jp
//


#include <iostream>
#include <vector>
using namespace std;


template <class Abel> struct BIT {
    vector<Abel> dat[2];
    Abel UNITY_SUM = 0;                     // to be set
    
    /* [1, n] */
    BIT(int n) { init(n); }
    void init(int n) { for (int iter = 0; iter < 2; ++iter) dat[iter].assign(n + 1, UNITY_SUM); }
    
    /* a, b are 1-indexed, [a, b) */
    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i += i & -i)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * -(a - 1)); sub_add(1, a, x); sub_add(0, b, x * (b - 1)); sub_add(1, b, x * (-1));
    }
    
    /* a is 1-indexed, [a, b) */
    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a; i > 0; i -= i & -i) res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b - 1) + sub_sum(1, b - 1) * (b - 1) - sub_sum(0, a - 1) - sub_sum(1, a - 1) * (a - 1);
    }
    
    /* debug */
    void print() {
        for (int i = 1; i < (int)dat[0].size(); ++i) cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};



int main() {
    int N, Q; cin >> N >> Q;
    BIT<long long> bit(N);
    for (int query = 0; query < Q; ++query) {
        int type; cin >> type;
        if (type == 0) {
            int s, t; long long x; cin >> s >> t >> x;
            bit.add(s, t+1, x);
        }
        else {
            int s, t; cin >> s >> t;
            cout << bit.sum(s, t+1) << endl;
        }
    }
}
