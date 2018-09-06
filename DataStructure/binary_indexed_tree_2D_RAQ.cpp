//
// Binary Indexed Tree 2D (「領域加算」「領域和取得」両対応)
//
// verified:
//   Codeforces 198 DIV1 D - Iahub and Xors
//     http://codeforces.com/problemset/problem/341/D
//

#include <iostream>
#include <vector>
#include <queue>
using namespace std;


template <class Abel> struct BIT2D {
    const Abel UNITY_SUM = 0;						// to be set
    vector<vector<Abel> > dat[2][2];
    
    BIT2D(int n = 1, int m = 1) { init(n, m); }
    void init(int n = 1, int m = 1) {
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) dat[i][j].resize(n+1, vector<Abel>(m+1, UNITY_SUM));
    }
    
    /* a1 <= x < b1, a2 <= y < b2, 1-indexd */
    inline void subsub_add(int f, int s, int a, int b, Abel x) {
        for (int i = a; i < (int)dat[f][s].size(); i += i & -i)
            for (int j = b; j < (int)dat[f][s][0].size(); j += j & -j)
                dat[f][s][i][j] = dat[f][s][i][j] + x;
    }
    inline void sub_add(int a1, int a2, Abel x) {
        subsub_add(0, 0, a1, a2, x*(a1-1)*(a2-1));
        subsub_add(1, 0, a1, a2, -x*(a1-1));
        subsub_add(0, 1, a1, a2, -x*(a2-1));
        subsub_add(1, 1, a1, a2, x);
    }
    inline void add(int a1, int a2, int b1, int b2, Abel x) {
        sub_add(a1, a2, x);
        sub_add(a1, b2, -x);
        sub_add(b1, a2, -x);
        sub_add(b1, b2, x);
    }
    
    /* a1 <= x < b1, a2 <= y < b2, 1-indexd */
    inline Abel subsub_sum(int f, int s, int a, int b) {
        Abel res = 0;
        for (int i = a; i > 0; i -= i & -i)
            for (int j = b; j > 0; j -= j & -j)
                res = res + dat[f][s][i][j];
        return res;
    }
    inline Abel sub_sum(int a1, int a2) {
        Abel res = 0;
        res += subsub_sum(0, 0, a1-1, a2-1);
        res += subsub_sum(1, 0, a1-1, a2-1) * (a2-1);
        res += subsub_sum(0, 1, a1-1, a2-1) * (a1-1);
        res += subsub_sum(1, 1, a1-1, a2-1) * (a1-1) * (a2-1);
        return res;
    }
    inline Abel sum(int a1, int a2, int b1, int b2) {
        return sub_sum(b1, b2) - sub_sum(a1, b2) - sub_sum(b1, a2) + sub_sum(a1, a2);
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



/* xor に対応するため */
struct XOR {
    long long val;
    
    XOR() : val(0) {}
    XOR(long long val_) { this->val = val_; }
    XOR operator = (long long val_) { this->val = val_; return *this; }
    inline XOR operator - () { return val; }
    inline const XOR& operator += (const XOR &x);
    inline const XOR& operator -= (const XOR &x);
};
inline XOR operator + (XOR x, XOR y) { return (x.val ^ y.val); }
inline XOR operator - (XOR x, XOR y) { return (x.val ^ y.val); }
inline XOR operator * (XOR x, int p) { if (p & 1) { return x; } else { return 0; } }
inline const XOR& XOR::operator += (const XOR &x) { *this = *this + x; return *this; }
inline const XOR& XOR::operator -= (const XOR &x) { *this = *this - x; return *this; }



int main() {
    int n, m;
    int q, x0, y0, x1, y1;
    long long v;
    
    cin >> n >> m;
    BIT2D<XOR> bit(n, n);
    for (int i = 0; i < m; ++i) {
        scanf("%d", &q);
        if (q == 2) {
            scanf("%d %d %d %d %lld", &x0, &y0, &x1, &y1, &v);
            bit.add(x0, y0, x1+1, y1+1, v);
        }
        if (q == 1) {
            scanf("%d %d %d %d", &x0, &y0, &x1, &y1);
            cout << bit.sum(x0, y0, x1+1, y1+1).val << endl;
        }
    }
}

