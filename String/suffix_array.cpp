//
// Suffix Array (SA-IS) (O(N))
//   G. Nong, S. Zhang, and W. H. Chan:
//   Two Efficient Algorithms for Linear Time Suffix Array Construction
//
// verified (suffix array の lcp を sparse table で求める)
//
//   Yosupo Judge Suffix Array
//     https://judge.yosupo.jp/problem/suffixarray
//
//   ARC 060 F - 最良表現
//     https://beta.atcoder.jp/contests/arc060/tasks/arc060_d
//
//   AOJ 2644 Longest Match
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2644
//


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
        int n = (int)vec.size(), h = 0;
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

// SA-IS (O(N))
template<class Str> struct SuffixArray {
    // data
    Str str;
    vector<int> sa;    // sa[i] : the starting index of the i-th smallest suffix (i = 0, 1, ..., n)
    vector<int> rank;  // rank[sa[i]] = i
    vector<int> lcp;   // lcp[i]: the lcp of sa[i] and sa[i+1] (i = 0, 1, ..., n-1)
    SparseTable<int> st;  // use for calcultating lcp(i, j)

    // getter
    int& operator [] (int i) {
        return sa[i];
    }
    vector<int> get_sa() { return sa; }
    vector<int> get_rank() { return rank; }
    vector<int> get_lcp() { return lcp; }

    // constructor
    SuffixArray(const Str& str_) : str(str_) {
        build_sa();
    }
    void init(const Str& str_) {
        str = str_;
        build_sa();
    }
    void build_sa() {
        vector<int> s;
        for (int i = 0; i < (int)str.size(); ++i) {
            s.push_back(str[i] + 1);
        }
        s.push_back(0);
        sa = sa_is(s);
        calcLCP(s);
        buildSparseTable();
    }

    // SA-IS
    // upper: # of characters 
    vector<int> sa_is(vector<int> &s, int upper = 256) {
        int N = (int)s.size();
        if (N == 0) return {};
        else if (N == 1) return {0};
        else if (N == 2) {
            if (s[0] < s[1]) return {0, 1};
            else return {1, 0};
        }

        vector<int> isa(N);
        vector<bool> ls(N, false);
        for (int i = N - 2; i >= 0; --i) {
            ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
        }
        vector<int> sum_l(upper + 1, 0), sum_s(upper + 1, 0);
        for (int i = 0; i < N; ++i) {
            if (!ls[i]) ++sum_s[s[i]];
            else ++sum_l[s[i] + 1];
        }
        for (int i = 0; i <= upper; ++i) {
            sum_s[i] += sum_l[i];
            if (i < upper) sum_l[i + 1] += sum_s[i];
        }

        auto induce = [&](const vector<int> &lms) -> void {
            fill(isa.begin(), isa.end(), -1);
            vector<int> buf(upper + 1);
            copy(sum_s.begin(), sum_s.end(), buf.begin());
            for (auto d: lms) {
                if (d == N) continue;
                isa[buf[s[d]]++] = d;
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            isa[buf[s[N - 1]]++] = N - 1;
            for (int i = 0; i < N; ++i) {
                int v = isa[i];
                if (v >= 1 && !ls[v - 1]) {
                    isa[buf[s[v - 1]]++] = v - 1;
                }
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            for (int i = N - 1; i >= 0; --i) {
                int v = isa[i];
                if (v >= 1 && ls[v - 1]) {
                    isa[--buf[s[v - 1] + 1]] = v - 1;
                }
            }
        };
            
        vector<int> lms, lms_map(N + 1, -1);
        int M = 0;
        for (int i = 1; i < N; ++i) {
            if (!ls[i - 1] && ls[i]) {
                lms_map[i] = M++;
            }
        }
        lms.reserve(M);
        for (int i = 1; i < N; ++i) {
            if (!ls[i - 1] && ls[i]) {
                lms.push_back(i);
            }
        }
        induce(lms);

        if (M) {
            vector<int> lms2;
            lms2.reserve(isa.size());
            for (auto v: isa) {
                if (lms_map[v] != -1) lms2.push_back(v);
            }
            int rec_upper = 0;
            vector<int> rec_s(M);
            rec_s[lms_map[lms2[0]]] = 0;
            for (int i = 1; i < M; ++i) {
                int l = lms2[i - 1], r = lms2[i];
                int nl = (lms_map[l] + 1 < M) ? lms[lms_map[l] + 1] : N;
                int nr = (lms_map[r] + 1 < M) ? lms[lms_map[r] + 1] : N;
                bool same = true;
                if (nl - l != nr - r) same = false;
                else {
                    while (l < nl) {
                        if (s[l] != s[r]) break;
                        ++l, ++r;
                    }
                    if (l == N || s[l] != s[r]) same = false;
                }
                if (!same) ++rec_upper;
                rec_s[lms_map[lms2[i]]] = rec_upper;
            }
            auto rec_sa = sa_is(rec_s, rec_upper);

            vector<int> sorted_lms(M);
            for (int i = 0; i < M; ++i) {
                sorted_lms[i] = lms[rec_sa[i]];
            }
            induce(sorted_lms);
        }
        return isa;
    }

    // find min id that str.substr(sa[id]) >= T
    int lower_bound(const Str& T) {
        int left = -1, right = sa.size();
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (str.compare(sa[mid], string::npos, T) < 0)
                left = mid;
            else
                right = mid;
        }
        return right;
    }

    // find min id that str.substr(sa[id], T.size()) > T
    int upper_bound(const Str& T) {
        int left = -1, right = sa.size();
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (str.compare(sa[mid], T.size(), T) <= 0)
                left = mid;
            else
                right = mid;
        }
        return right;
    }

    // search
    bool is_contain(const Str& T) {
        int lb = lower_bound(T);
        if (lb >= sa.size()) return false;
        return str.compare(sa[lb], T.size(), T) == 0;
    }

    // find lcp
    void calcLCP(const vector<int> &s) {
        int N = (int)s.size();
        rank.assign(N, 0), lcp.assign(N, 0);
        for (int i = 0; i < N; ++i) rank[sa[i]] = i;
        int h = 0;
        for (int i = 0; i < N - 1; ++i) {
            int pi = sa[rank[i] - 1];
            if (h > 0) --h;
            for (; pi + h < N && i + h < N; ++h) {
                if (s[pi + h] != s[i + h]) break;
            }
            lcp[rank[i] - 1] = h;
        }
    }
    
    // build sparse table for calculating lcp
    void buildSparseTable() {
        st.init(lcp);
    }

    // calc lcp of str.sutstr(a) and str.substr(b)
    int getLCP(int a, int b) {
        return st.get(min(rank[a], rank[b]), max(rank[a], rank[b]));
    }
};


//////////////////////////////////////
// Examples
//////////////////////////////////////

// ARC 060 F
void solveARC060F() {
    string str; 
    cin >> str;
    int n = (int)str.size();

    // n の約数を求める
    vector<long long> divs;
    for (long long i = 1LL; i*i <= n; ++i) {
        if (n%i == 0LL) {
            divs.push_back(i);
            long long temp = n/i;
            if (i != temp) divs.push_back(temp);
        }
    }
    sort(divs.begin(), divs.end());

    // 求める
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

// AOJ 2644
void solveAOJ2644() {
    // Suffix Array の構築
    string S;
    cin >> S;
    SuffixArray<string> suf(S);
    vector<int> sa = suf.get_sa();

    // Suffix Array の区間に関する RMQ
    SparseTable<int> min_st(sa);
    for (auto& v : sa) v = -v;
    SparseTable<int> max_st(sa);

    auto solve = [&](const string& x, const string& y) -> int {
        if (!suf.is_contain(x) || !suf.is_contain(y)) {
            return 0;
        }
        int xl = suf.lower_bound(x);
        int xr = suf.upper_bound(x);
        int yl = suf.lower_bound(y);
        int yr = suf.upper_bound(y);
        int xres = min_st.get(xl, xr);
        int yres = -max_st.get(yl, yr);

        if (xres > yres || xres + x.size() > yres + y.size()) 
            return 0;
        else 
            return yres - xres + (int)y.size();   
    };

    // クエリ処理
    int Q;
    cin >> Q;
    while (Q--) {
        string x, y;
        cin >> x >> y;
        cout << solve(x, y) << endl;
    }
}

int main() {
    solveARC060F();
    //solveAOJ2644();
}