//
// Manacher のアルゴリズム
//   res[i] := S[i] を中心とする S 上の奇数長の最長回文の半径 (中心を含まない)
//   偶数長の回文も求めるために S = "abc" に対して "a$b$c" といった処理をする
//
// verified:
//   RUPC 2019 day3-E 往復文字列
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2934
//
//   ウクーニャたんお誕生日コンテスト D - ukuku
//     https://atcoder.jp/contests/ukuku09/tasks/ukuku09_d
//


#include <bits/stdc++.h>
using namespace std;


// Manacher
struct Manacher {
    string S;
    vector<int> radius_odd, radius_even;

    // construct
    Manacher(const string &S_) : S(S_) { init(S); }
    void init(const string &S_) {
        S = S_;
        string S2 = "";
        for (int i = 0; i < (int)S.size(); ++i) {
            S2 += S[i];
            if (i+1 < (int)S.size()) S2 += "$";
        }
        construct(S2);
    }
    vector<int> construct(const string &S2) {
        vector<int> len(S2.size());
        int i = 0, j = 0;
        while (i < (int)S2.size()) {
            while (i-j >= 0 && i+j < (int)S2.size() && S2[i-j] == S2[i+j]) ++j;
            len[i] = j;
            int k = 1;
            while (i-k >= 0 && i+k < (int)S2.size() && k+len[i-k] < j) {
                len[i+k] = len[i-k];
                ++k;
            }
            i += k, j -= k;
        }
        radius_odd.assign(S.size(), 0), radius_even.assign(S.size()+1, 0);
        for (int i = 0; i < (int)S.size(); ++i) {
            radius_odd[i] = (len[i*2] - 1) / 2;
            if (i > 0) radius_even[i] = len[i*2-1] / 2;
        }
        return len;
    }

    // radius, center is i (0 <= i < N)
    int get_odd(int i) { return radius_odd[i]; }

    // radius, center is between i-1 and i (0 <= i <= N)
    int get_even(int i) { return radius_even[i]; }

    // judge if [left, right) is palindrome
    bool is_palindrome(int left, int right) {
        int mid = (left + right) / 2;
        if ((right - left) & 1) return ( get_odd(mid) == (right - left + 1)/2);
        else return (get_even(mid) == (right - left)/2);
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void RUPC_2019_day3_E() {
    int N;
    string S;
    cin >> N >> S;

    Manacher m(S);
    int res = N;
    for (int len = 2; len < N; ++len) {
        bool ok = true;
        for (int i = len - 1; i < N; i += len - 1) {
            if (m.get_odd(i) != min(i, N-i-1)) ok = false;
        }
        if (ok) {
            res = len;
            break;
        }
    }
    cout << res << endl;
}

// Sparse Table
template<class MeetSemiLattice> struct SparseTable {
    vector<vector<MeetSemiLattice> > dat;
    vector<int> height;
    
    SparseTable() { }
    SparseTable(const vector<MeetSemiLattice> &vec) { init(vec); }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) < n) ++h;
        dat.assign(h, vector<MeetSemiLattice>(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i)
            for (int j = 0; j < n; ++j)
                dat[i][j] = max(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
    }
    
    MeetSemiLattice get(int a, int b) {
        if (a >= b) return -1;
        return max(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};

void UKUKU_D() {
    int N, Q;
    string S;
    cin >> N >> Q >> S;
    Manacher mana(S);
    SparseTable<int> st(mana.radius_odd);
    while (Q--) {
        int l, r;
        cin >> l >> r;
        --l;
        
        // 半径 r 以下の回文が存在するような最大の r を求める
        int low = 0, high = (r-l)/2+1;
        while (high - low > 1) {
            int mid = (low + high) / 2;
            if (st.get(l+mid, r-mid) >= mid) low = mid;
            else high = mid;
        }
        cout << low * 2 + 1 << endl;
    }
}


int main() {
    //RUPC_2019_day3_E();
    UKUKU_D();
}

