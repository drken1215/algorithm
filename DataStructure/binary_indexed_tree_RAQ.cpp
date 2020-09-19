//
// Binary Indexed Tree (「区間加算」「区間和取得」両対応)
//
// verified:
//   AOJ Course DSL_2_G RSQ and RAQ
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_G&lang=ja
//


#include <iostream>
#include <vector>
using namespace std;


template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat[2];

    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity) {
        init(n);
    }
    void init(int n) {
        for (int iter = 0; iter < 2; ++iter)
            dat[iter].assign(n + 1, UNITY_SUM);
    }
    
    // [a, b), a and b are 0-indexed
    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i |= i + 1)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * (-a));
        sub_add(1, a, x);
        sub_add(0, b, x * b);
        sub_add(1, b, x * (-1));
    }
    
    // [a, b), a and b are 0-indexed
    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b)
            + sub_sum(1, b) * b
            - sub_sum(0, a)
            - sub_sum(1, a) * a;
    }
    
    // debug
    void print() {
        for (int i = 0; i < (int)dat[0].size(); ++i)
            cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};



int main() {
    int N, Q;
    cin >> N >> Q;
    BIT<long long> bit(N);
    for (int query = 0; query < Q; ++query) {
        int type, s, t;
        cin >> type >> s >> t;
        --s;
        if (type == 0) {
            long long x;
            cin >> x;
            bit.add(s, t, x);
        }
        else {
            cout << bit.sum(s, t) << endl;
        }
    }
}
