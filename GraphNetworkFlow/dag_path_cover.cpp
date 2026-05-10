//
// DAG の最小パス被覆
//
// verified
//   ABC 237 Ex - Hakata
//     https://atcoder.jp/contests/abc237/tasks/abc237_h
//


#include <bits/stdc++.h>
using namespace std;


// DAG min path-cover by Hopcroft-Karp
struct DagPathCover {
    const int NOT_MATCHED = -1;
    
    // input
    int V;
    vector<vector<int>> list; // left to right

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<vector<int>> rlist; // right to left
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    DagPathCover(int V) : V(V), list(V), rlist(V) {}
    void add_edge(int from, int to) {
        assert(from >= 0 && from < V);
        assert(to >= 0 && to < V);
        list[from].emplace_back(to);
        rlist[to].emplace_back(from);
    }

    // getter, debugger
    vector<int> &operator [] (int i) { return list[i]; }
    const vector<int> &operator [] (int i) const { return list[i]; }
    constexpr int size() const { return V; }
    friend ostream& operator << (ostream& s, const DagPathCover& G) {
        s << endl;
        for (int i = 0; i < (int)G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < (int)G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != (int)G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // solver
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < V; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[V] = V;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int right : list[left]) {
                int next = rl[right];
                if (next == NOT_MATCHED) next = V;
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == V) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int right : list[left]) {
            int next = rl[right];
            if (next == NOT_MATCHED) next = V;
            if (level[next] > level[left] && hodfs(next)) {
                lr[left] = right;
                rl[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(V, false);
        matched.assign(V, false);
        level.assign(V + 1, -1);
        lr.assign(V, -1);
        rl.assign(V, -1);
        int max_matching = 0;
        while (true) {
            hobfs();
            seen.assign(V, false);
            bool finished = true;
            for (int left = 0; left < V; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++max_matching;
                    finished = false;
                }
            }
            if (finished) break;
        }
        for (int r = 0; r < V; r++) {
            if (rl[r] != NOT_MATCHED) lr[rl[r]] = r;
        }
        return V - max_matching;
    }

    // max stable set
    vector<int> get_stable_set() {
        vector<int> res;
        for (int v = 0; v < V; v++) {
            if (rl[v] == NOT_MATCHED) res.emplace_back(v);
        }
        return res;
    }

    // min path cover
    vector<vector<int>> get_path_cover() {
        auto srcs = get_stable_set();
        vector<vector<int>> res;
        for (auto s : srcs) {
            vector<int> path;
            int v = s;
            while (v != NOT_MATCHED) {
                path.emplace_back(v);
                v = lr[v];
            }
            res.emplace_back(path);
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 237 Ex - Hakata
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP = [](const MeetSemiLattice &l, const MeetSemiLattice &r) {
        return min(l, r);
    };
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() {}
    SparseTable(const vector<MeetSemiLattice> &vec) {
        init(vec);
    }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  {
        init(vec, op);
    }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<MeetSemiLattice>(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            for (int j = 0; j < n; ++j)
                dat[i][j] = OP(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
        }
    }
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
        init(vec);
    }
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};
template<class Str = string> struct SuffixArray {
    // data
    Str str;
    vector<int> sa;    // sa[i] : the starting index of the i-th smallest suffix (i = 0, 1, ..., n)
    vector<int> rank;  // rank[sa[i]] = i
    vector<int> lcp;   // lcp[i]: the lcp of sa[i] and sa[i+1] (i = 0, 1, ..., n-1)
    SparseTable<int> st;  // use for calcultating lcp(i, j)

    // getter
    int& operator [] (int i) { return sa[i]; }
    const int& operator [] (int i) const { return sa[i]; }
    vector<int> get_sa() { return sa; }
    vector<int> get_rank() { return rank; }
    vector<int> get_lcp() { return lcp; }

    // constructor
    SuffixArray() {}
    SuffixArray(const Str& str_, bool no_limit_elements = false) : str(str_) {
        build_sa(no_limit_elements);
    }
    void init(const Str& str_, bool no_limit_elements = false) {
        str = str_;
        build_sa(no_limit_elements);
    }
    void build_sa(bool no_limit_elements = false) {
        vector<int> s;
        int num_of_chars = 256;
        if (!no_limit_elements) {
            for (int i = 0; i < (int)str.size(); ++i) {
                s.push_back(str[i] + 1);
            }
        } else {
            unordered_map<int,int> dict;
            for (int i = 0; i < (int)str.size(); ++i) {
                if (!dict.count(str[i])) dict[str[i]] = dict.size();
            }
            for (int i = 0; i < (int)str.size(); ++i) {
                s.push_back(dict[str[i]] + 1);
            }
            num_of_chars = (int)dict.size();
        }
        s.push_back(0);
        sa = sa_is(s, num_of_chars);
        build_lcp(s);
        build_sparse_table();
    }

    // SA-IS
    // num_of_chars: # of characters
    vector<int> sa_is(vector<int> &s, int num_of_chars) {
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
        vector<int> sum_l(num_of_chars + 1, 0), sum_s(num_of_chars + 1, 0);
        for (int i = 0; i < N; ++i) {
            if (!ls[i]) ++sum_s[s[i]];
            else ++sum_l[s[i] + 1];
        }
        for (int i = 0; i <= num_of_chars; ++i) {
            sum_s[i] += sum_l[i];
            if (i < num_of_chars) sum_l[i + 1] += sum_s[i];
        }

        auto induce = [&](const vector<int> &lms) -> void {
            fill(isa.begin(), isa.end(), -1);
            vector<int> buf(num_of_chars + 1);
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
    
    // find min id that sa[id] >= str.substr(l, r-l)
    int lower_bound(int l, int r) {
        int left = -1, right = rank[l];
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (st.get(mid, rank[l]) < r - l) left = mid;
            else right = mid;
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
    void build_lcp(const vector<int> &s) {
        int N = (int)s.size();
        rank.assign(N, 0), lcp.assign(N - 1, 0);
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
    void build_sparse_table() {
        st.init(lcp);
    }

    // calc lcp of str.sutstr(a) and str.substr(b)
    int get_lcp(int a, int b) {
        return st.get(min(rank[a], rank[b]), max(rank[a], rank[b]));
    }

    // debug
    void dump() {
        for (int i = 0; i < sa.size(); ++i) {
            cout << i << ": " << sa[i] << ", " << str.substr(sa[i]) << endl;
        }
    }
};
void ABC_237_Ex() {
    string S;
    cin >> S;
    int N = S.size();

    set<string> ins;
    for (int l = 0; l < N; l++) for (int r = l+1; r <= N; r++) {
        bool ok = true;
        for (int i = l, j = r-1; i < j; i++, j--) if (S[i] != S[j]) ok = false;
        if (ok) ins.insert(S.substr(l, r-l));
    }
    vector<string> alts;
    for (auto s : ins) alts.emplace_back(s);
    int V = alts.size();

    string T = "";
    vector<int> starts(V, 0);
    for (int v = 0; v < V; v++) {
        starts[v] = (int)T.size();
        T += alts[v] + "$";
    }
    SuffixArray<string> sa(T);

    // does i include j ?
    auto check = [&](int i, int j) -> bool {
        if (alts[i].size() <= alts[j].size()) return false;
        int diff = (int)alts[i].size() - (int)alts[j].size();
        for (int k = 0; k <= diff; k++) {
            if (sa.get_lcp(starts[i]+k, starts[j]) >= (int)alts[j].size()) 
                return true;
        }
        return false;
    };

    DagPathCover G(V);
    for (int v = 0; v < V; v++) {
        for (int v2 = 0; v2 < V; v2++) {
            if (v == v2) continue;
            if (check(v, v2)) G.add_edge(v, v2);
        }
    }
    int res = G.solve();
    cout << res << endl;
}


int main() {
    ABC_237_Ex();
}