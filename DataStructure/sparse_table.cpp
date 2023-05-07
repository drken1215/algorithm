//
// RMQ (by sparse table)
//
// verified (suffix array の lcp を sparse table で求める):
//   ARC 060 F - 最良表現
//     https://beta.atcoder.jp/contests/arc060/tasks/arc060_d
//


/*
  動的更新はできないが、最初に静的に構築してしまえば、範囲最小値を O(1) でとり出せる
    交叉半束 (meet-semi-lattice) 上で動作
    半群 (semi-group) 上で動作する上位互換バージョン (disjoint sparse table) もある

  init(vec): 配列を vec で初期化構築, O(n logn)
  get(a, b): [a, b) の最小値を取得
*/


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;


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
                dat[i][j] = min(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
    }
    
    MeetSemiLattice get(int a, int b) {
        return min(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};


// Suffix Array ( Manber&Myers: O(n (logn)^2) )
struct SuffixArray {
    string str;
    vector<int> sa;         // sa[i] : the starting index of the i-th smallest suffix (i = 0, 1, ..., n)
    vector<int> lcp;        // lcp[i]: the lcp of sa[i] and sa[i+1] (i = 0, 1, ..., n-1)
    inline int& operator [] (int i) {return sa[i];}
    
    SuffixArray(const string& str_) : str(str_) { buildSA(); calcLCP(); }
    void init(const string& str_) { str = str_; buildSA(); calcLCP(); }
    
    // build SA
    vector<int> rank_sa, tmp_rank_sa;
    struct CompareSA {
        int n, k;
        const vector<int> &rank;
        CompareSA(int n, int k, const vector<int> &rank_sa) : n(n), k(k), rank(rank_sa) {}
        bool operator()(int i, int j) {
            if (rank[i] != rank[j]) return (rank[i] < rank[j]);
            else {
                int rank_ik = (i + k <= n ? rank[i + k] : -1);
                int rank_jk = (j + k <= n ? rank[j + k] : -1);
                return (rank_ik < rank_jk);
            }
        }
    };
    void buildSA() {
        int n = (int)str.size();
        sa.resize(n+1), lcp.resize(n+1), rank_sa.resize(n+1), tmp_rank_sa.resize(n+1);
        for (int i = 0; i < n; ++i) sa[i] = i, rank_sa[i] = (int)str[i];
        sa[n] = n, rank_sa[n] = -1;
        for (int k = 1; k <= n; k *= 2) {
            CompareSA csa(n, k, rank_sa);
            sort(sa.begin(), sa.end(), csa);
            tmp_rank_sa[sa[0]] = 0;
            for (int i = 1; i <= n; ++i) {
                tmp_rank_sa[sa[i]] = tmp_rank_sa[sa[i - 1]];
                if (csa(sa[i - 1], sa[i])) ++tmp_rank_sa[sa[i]];
            }
            for (int i = 0; i <= n; ++i) rank_sa[i] = tmp_rank_sa[i];
        }
    }
    vector<int> rsa;
    SparseTable<int> st;
    void calcLCP() {
        int n = (int)str.size();
        rsa.resize(n+1);
        for (int i = 0; i <= n; ++i) rsa[sa[i]] = i;
        lcp.resize(n+1);
        lcp[0] = 0;
        int cur = 0;
        for (int i = 0; i < n; ++i) {
            int pi = sa[rsa[i] - 1];
            if (cur > 0) --cur;
            for (; pi + cur < n && i + cur < n; ++cur) {
                if (str[pi + cur] != str[i + cur]) break;
            }
            lcp[rsa[i] - 1] = cur;
        }
        st.init(lcp);
    }
    
    // calc lcp
    int getLCP(int a, int b) {          // lcp of str.sutstr(a) and str.substr(b)
        return st.get(min(rsa[a], rsa[b]), max(rsa[a], rsa[b]));
    }
};


vector<long long> divisor(long long n) {
    vector<long long> res;
    for (long long i = 1LL; i*i <= n; ++i) {
        if (n%i == 0LL) {
            res.push_back(i);
            long long temp = n/i;
            if (i != temp) res.push_back(temp);
        }
    }
    sort(res.begin(), res.end());
    return res;
}



int main() {
    string str; cin >> str;
    int n = (int)str.size();
    vector<long long> divs = divisor(n);
    long long syuuki = n;
    for (auto d : divs) {
        bool ok = true;
        for (int j = 0; j + d < n; ++j) {
            if (str[j] != str[j+d]) ok = false;
        }
        if (ok) syuuki = min(syuuki, d);
    }
    if (syuuki == n) cout << 1 << endl << 1 << endl;
    else if (syuuki == 1) cout << n << endl << 1 << endl;
    else {
        vector<int> cannot_cut(n*2, 0);
        SuffixArray sa(str);
        for (int d = 1; d < n; ++d) {
            for (int dd = d; dd < n; dd += d) {
                int lcp = sa.getLCP(0, dd);
                if (lcp < d) break;
                cannot_cut[dd + d] = true;
            }
            for (int dd = n-d*2; dd >= 0; dd -= d) {
                int lcp = sa.getLCP(dd, n-d);
                if (lcp < d) break;
                cannot_cut[dd] = true;
            }
        }
        int con = 0;
        for (int i = 1; i < n; ++i) if (!cannot_cut[i]) ++con;
        cout << 2 << endl << con << endl;
    }
}
