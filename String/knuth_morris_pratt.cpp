//
// Knuth-Morris-Pratt
//   fail[i] := pat[0:i) の suffix と pat の prefix がどれだけ一致するか (長さ i 未満で)
//   O(N) で構築できる
//   これを利用して以下を求められたる
//     ・pat[0:i) の最小周期 ("abcabcab" は "abc" の繰り返しを区切ることでできるので 3)
//     ・文字列 S において S[i:] の prefix が pat と一致するような i をすべて求める
//
// verified:
//   ABC 150 F - Xor Shift
//     https://atcoder.jp/contests/abc150/tasks/abc150_f
//
//   ARC 060 F - 最良表現
//     https://atcoder.jp/contests/arc060/tasks/arc060_d
//


#include <iostream>
#include <vector>
using namespace std;


// T = string or vector<long long>
template<class T> struct KMP {
    T pat;
    vector<int> fail;

    // construct
    KMP(const T &p) { init(p); }
    void init(const T &p) {
        pat = p;
        int m = (int)pat.size();
        fail.assign(m+1, -1);
        for (int i = 0, j = -1; i < m; ++i) {
            while (j >= 0 && pat[i] != pat[j]) j = fail[j];
            fail[i+1] = ++j;
        }
    }

    // the period of S[0:i]
    int period(int i) { return i - fail[i]; }
    
    // the index i such that S[i:] has the exact prefix p
    vector<int> match(const T &S) {
        int n = (int)S.size(), m = (int)pat.size();
        vector<int> res;
        for (int i = 0, k = 0; i < n; ++i) {
            while (k >= 0 && S[i] != pat[k]) k = fail[k];
            ++k;
            if (k >= m) res.push_back(i - m + 1), k = fail[k];
        }
        return res;
    }
};



//------------------------------//
// Examples
//------------------------------//

void solve(int N, const vector<long long> &a, const vector<long long> &b) {
    vector<long long> da(N*2), db(N);
    for (int i = 0; i < N*2; ++i) {
        int j = i % N;
        int k = (i + 1) % N;
        da[i] = a[j] ^ a[k];
        if (i < N) db[i] = b[j] ^ b[k];
    }
    KMP<vector<long long> > kmp(db);
    auto v = kmp.match(da);
    for (auto k : v) {
        if (k < N) cout << k << " " << (a[k] ^ b[0]) << endl;
    }
}

int main() {
    int N;
    cin >> N;
    vector<long long> a(N), b(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    for (int i = 0; i < N; ++i) cin >> b[i];
    solve(N, a, b);
}
