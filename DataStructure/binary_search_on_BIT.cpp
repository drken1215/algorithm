//
// BIT 上二分探索 (以下のクエリに答えられる、ただし v <= N として O(N) のメモリが必要、時間計算量は O(log N))
//   集合 S に要素 v を追加する
//   集合 S から要素 v を除く
//   集合 S の k 番目に大きい要素を答える
//
// verified:
//   ARC 033 C - データ構造
//     https://beta.atcoder.jp/contests/arc033/tasks/arc033_3
//


#include <iostream>
#include <vector>
using namespace std;


// get(k): binary search on BIT
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

    /* k-th number (k is 0-indexed) */
    int get(long long k) {
        ++k;
        int res = 0;
        int N = 1; while (N < (int)dat.size()) N *= 2;
        for (int i = N / 2; i > 0; i /= 2) {
            if (res + i < (int)dat.size() && dat[res + i] < k) {
                k = k - dat[res + i];
                res = res + i;
            }
        }
        return res + 1;
    }
    
    /* debug */
    void print() {
        for (int i = 1; i < (int)dat.size(); ++i) cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};



int main() {
    BIT<int> bit(200000);
    int Q;
    cin >> Q;
    for (int i = 0; i < Q; ++i) {
        int T, X;
        cin >> T >> X;
        if (T == 1) bit.add(X, 1); // insert X
        else {
            int val = bit.get(X - 1); // X is 1-indexed
            cout << val << endl;
            bit.add(val, -1); // delete val
        }
    }
}
