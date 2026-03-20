//
// Inversion Number of array a
//   a_i > a_j となる (i < j) の組の個数を数えます
//
// verified
//   AOJ Course ALDS1_5_D Recursion / Divide and Conquer - The Number of Inversions
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_5_D&lang=jp
//


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


template <class Abel> struct BIT {
    const Abel UNITY_SUM = 0;                       // to be set
    vector<Abel> dat;
    
    /* [1, n] */
    BIT(int n) : dat(n + 1, UNITY_SUM) { }
    void init(int n) { dat.assign(n + 1, UNITY_SUM); }
    
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

template<class T> long long inversion_number(const vector<T> &a) {
    int n = (int)a.size();
    vector<T> b = a;
    sort(b.begin(), b.end());
    b.erase(unique(b.begin(), b.end()), b.end());
    
    long long res = 0;
    BIT<int> bit(n);
    for (int i = 0; i < n; ++i) {
        int order = lower_bound(b.begin(), b.end(), a[i]) - b.begin();
        ++order;
        res += bit.sum(order+1, n+1);
        bit.add(order, 1);
    }
    return res;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int n; cin >> n;
    vector<int> a(n); for (int i = 0; i < n; ++i) cin >> a[i];
    cout << inversion_number(a) << endl;
}
