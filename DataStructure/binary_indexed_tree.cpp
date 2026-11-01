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


// BIT
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat;
    
    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity), dat(n, unity) { }
    void init(int n) {
        dat.assign(n, UNITY_SUM);
    }
    int size() const {
        return (int)dat.size();
    }
    
    // a is 0-indexed
    inline void add(int a, Abel x) {
        for (int i = a; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, a), a is 0-indexed, [a, b), a and b are 0-indexed
    inline Abel sum(int a) const {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    inline Abel sum(int a, int b) const {
        return sum(b) - sum(a);
    }
    inline Abel operator [] (int i) const {
        return sum(i, i + 1);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const BIT &bit) {
        for (int i = 0; i < (int)bit.size(); ++i) s << bit[i] << " ";
        return s;
    }
};


//------------------------------//
// Examples
//------------------------------//

int main() {
    int N, Q; cin >> N >> Q;
    BIT<int> bit(N);
    for (int query = 0; query < Q; ++query) {
        int com, x, y;
        cin >> com >> x >> y;
        --x;
        if (com == 0) bit.add(x, y);
        else cout << bit.sum(x, y) << endl;
    }
}
