//
// Binary Trie
//
// verified:
//   Yosupo Library Checker - Set Xor-Min (for insert(), erase(), get_max(add))
//     https://judge.yosupo.jp/problem/set_xor_min
//
//   AtCoder ARC 033 C - データ構造 (for insert(), get_kth())
//     https://atcoder.jp/contests/arc033/tasks/arc033_3
//
//   AtCoder 第5回 ドワンゴからの挑戦状 本選 B - XOR Spread (for erase(), add_all(), get_min())
//     https://atcoder.jp/contests/dwacon5th-final/tasks/dwacon5th_final_b
//
//   ABC 281 E - Least Elements (for insert(), get_min(), get_max())
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
//
//   AOJ 3333 - Range Xor Query (for get_min(add))
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=3333
//
//   ABC 451 G - Minimum XOR Walk (for insert(), count_lower(add))
//     https://atcoder.jp/contests/abc451/tasks/abc451_g
//


#include <bits/stdc++.h>
using namespace std;


// Binary Trie
template<typename INT, size_t MAX_DIGIT> struct BinaryTrie {
    // inner values
    struct Node {
        size_t count;
        Node *prev, *left, *right;
        Node(Node *prev) : count(0), prev(prev), left(nullptr), right(nullptr) {}
    };
    INT lazy;
    Node *root;

    // constructor and getters
    BinaryTrie() : lazy(0), root(emplace(nullptr)) {}
    INT get(Node *v) {
        if (!v) return -1;
        INT res = 0;
        for (int i = 0; i < MAX_DIGIT; ++i) {
            if (!v) break;
            if (v == v->prev->right) res |= INT(1)<<i;
            v = v->prev;
        }
        return res ^ lazy;
    }
    Node* find(INT val) {
        INT nval = val ^ lazy;
        Node *v = root;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            bool flag = (nval >> i) & 1;
            if (flag) v = v->right;
            else v = v->left;
            if (!v) return v;
        }
        return v;
    }
    void add_all(INT add) { lazy ^= add; }

    // insert
    Node* emplace(Node *prev) {
        return new Node(prev);
    }
    void insert(INT val, size_t k = 1) {
        INT nval = val ^ lazy;
        Node *v = root;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            bool flag = (nval >> i) & 1;
            if (flag && !v->right) v->right = emplace(v);
            if (!flag && !v->left) v->left = emplace(v);
            if (flag) v = v->right;
            else v = v->left;
        }
        v->count += k;
        while ((v = v->prev)) v->count = get_count(v->left) + get_count(v->right);
    }
    
    // erase
    Node* clear(Node *v) {
        if (!v || get_count(v)) return v;
        delete(v);
        return nullptr;
    }
    bool erase(Node *v, size_t k = 1) {
        if (!v) return false;
        v->count -= k;
        while ((v = v->prev)) {
            v->left = clear(v->left);
            v->right = clear(v->right);
            v->count = get_count(v->left) + get_count(v->right);
        }
        return true;
    }
    bool erase(INT val, size_t k = 1) {
        auto v = find(val);
        return erase(v, k);
    }

    // counter
    size_t size() const { return get_count(root); }
    size_t get_count(Node *v) const { return v ? v->count : 0; }
    size_t count_lower(INT val, INT add = 0) {  // calc #{x | (x ^ add) < val}
        if (!root) return 0;
        INT addlazy = add ^ lazy;
        Node *v = root;
        size_t res = 0;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            Node *left = v->left, *right = v->right;
            if ((addlazy >> i) & 1) swap(left, right);
            bool flag = (val >> i) & 1;
            if (flag) {
                res += get_count(left);
                v = right;
            }
            else v = left;
        }
        return res;
    }
    size_t count(INT left, INT right) {
        return count_lower(right) - count_lower(left);
    }
    size_t count(INT val) {
        return count(val, val + 1);
    }

    // max (with xor-addition of val) and min (with xor-addition of add)
    Node* get_max_node(INT add = 0) {
        INT addlazy = add ^ lazy;
        Node* v = root;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            bool flag = (addlazy >> i) & 1;
            if (!v->right) v = v->left;
            else if (!v->left) v = v->right;
            else if (flag) v = v->left;
            else v = v->right;
        }
        return v;
    }
    INT get_max(INT add = 0) {
        return get(get_max_node(add)) ^ add;
    }
    Node* get_min_node(INT add = 0) {
        return get_max_node(~add & ((INT(1)<<MAX_DIGIT)-1));
    }
    INT get_min(INT add = 0) {
        return get(get_min_node(add)) ^ add;
    }

    // find k-th val, k is 0-indexed
    Node* get_kth_node(size_t k, INT add = 0) {
        INT addlazy = add ^ lazy;
        Node *v = root;
        if (get_count(v) <= k) return nullptr;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            bool flag = (addlazy >> i) & 1;
            Node *left = (flag ? v->right : v->left);
            Node *right = (flag ? v->left : v->right);
            if (get_count(left) <= k) k -= get_count(left), v = right;
            else v = left;
        }
        return v;
    }
    INT get_kth(size_t k, INT add = 0) {
        return get(get_kth_node(k, add));
    }

    // next (including val)
    INT next(INT val) {
        return get(lower_bound(val));
    }
   
    // lower_bound (with xor-addition of add), upper_bound (with xor-addition of add)
    Node* get_cur_node(Node *v, int i) {
        if (!v) return v;
        Node *left = v->left, *right = v->right;
        if ((lazy >> i) & 1) swap(left, right);
        if (left) return get_cur_node(left, i+1);
        else if (right) return get_cur_node(right, i+1);
        return v;
    }
    Node* get_next_node(Node *v, int i) {
        if (!v->prev) return nullptr;
        Node *left = v->prev->left, *right = v->prev->right;
        if ((lazy >> (i+1)) & 1) swap(left, right);
        if (v == left && right) return get_cur_node(right, i);
        else return get_next_node(v->prev, i+1);
    }
    Node* lower_bound(INT val, INT add = 0) {
        INT addlazy = add ^ lazy;
        INT nval = val ^ addlazy;
        Node *v = root;
        for (int i = MAX_DIGIT-1; i >= 0; --i) {
            if (!v) break;
            bool flag = (nval >> i) & 1;
            if (flag && v->right) v = v->right;
            else if (!flag && v->left) v = v->left;
            else if ((val >> i) & 1) return get_next_node(v, i);
            else return get_cur_node(v, i);
        }
        return v;
    }
    Node* upper_bound(INT val, INT add = 0) {
        return lower_bound(val + 1, add);
    }

    // debug
    void print(Node *v, string prefix = "") {
        if (!v) return;
        cout << prefix << ": " << v->count << endl;
        print(v->left, prefix + "0");
        print(v->right, prefix + "1");
    }
    void print() {
        print(root);
    }
    vector<INT> eval(Node *v, int digit) const {
        vector<INT> res;
        if (!v) return res;
        if (!v->left && !v->right) {
            for (int i = 0; i < get_count(v); ++i) res.push_back(0);
            return res;
        }
        const auto& left = eval(v->left, digit-1);
        const auto& right = eval(v->right, digit-1);
        for (auto val : left) res.push_back(val);
        for (auto val : right) res.push_back(val + (INT(1)<<digit));
        return res;
    }
    vector<INT> eval() const {
        auto res = eval(root, MAX_DIGIT-1);
        for (auto &val : res) val ^= lazy;
        return res;
    }
    friend ostream& operator << (ostream &os,
                                 const BinaryTrie<INT, MAX_DIGIT> &bt) {
        auto res = bt.eval();
        for (auto val : res) os << val << " ";
        return os;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Set Xor-Min
void YosupoSetXorMin() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int Q;
    cin >> Q;
    BinaryTrie<int, 30> bt;
    while (Q--) {
        int type, x;
        cin >> type >> x;
        if (type == 0) {
            if (!bt.count(x)) bt.insert(x);
        } else if (type == 1) {
            if (bt.count(x)) bt.erase(x);
        } else {
            cout << bt.get_min(x) << endl;
        }
    }
}

// AtCoder ARC 033 C - データ構造
void ARC033_C() {
    int Q;
    cin >> Q;
    BinaryTrie<int, 30> bt;
    while (Q--) {
        int type, x;
        cin >> type >> x;
        if (type == 1) bt.insert(x);
        else {
            --x;
            auto v = bt.get_kth(x);
            cout << v << endl;
            bt.erase(v);
        }
    }
}

// AtCoder 第5回 ドワンゴからの挑戦状 本選 B - XOR Spread
void AtCoderDwango5Honsen_B() {
    int N;
    cin >> N;
    vector<int> a(N), s(N), res(N, 0);
    for (int i = 0; i < N; ++i) {
        cin >> a[i];
        if (i == 0) s[i] = a[0];
        else s[i] = s[i-1] ^ a[i];
    }
    BinaryTrie<int, 30> bt;
    for (int i = 0; i < N-1; ++i) bt.insert(s[i]);
    for (int i = 0; i < N-1; ++i) {
        res[i] = bt.get_min();
        bt.erase(res[i]);
        bt.add_all(res[i]);
    }
    res.back() = s.back() ^ bt.lazy;
    for (int i = 0; i < N; ++i) cout << res[i] << " ";
    cout << endl;
}

// ABC 281 E - Least Elements
void ABC281_E() {
    // 入力
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 小さい順に K 個の総和
    long long sum = 0;

    // Binary Trie
    BinaryTrie<int, 30> left, right;

    // push と pop
    auto push = [&](long long x) -> void {
        // とりあえず x を left に挿入する
        left.insert(x);
        sum += x;

        // left のサイズが K を超えるなら left の最大値を right に移す
        if (left.size() > K) {
            long long y = left.get_max();
            sum -= y;
            left.erase(y);
            right.insert(y);
        }
    };
    auto pop = [&](long long x) -> void {
        // とりあえず x を削除する
        if (x <= left.get_max()) {
            left.erase(x);
            sum -= x;
        } else {
            right.erase(x);
        }

        // left のサイズが K 未満になるなら right の最小値を left に移す
        if (left.size() < K) {
            long long y = right.get_min();
            sum += y;
            left.insert(y);
            right.erase(y);
        }
    };

    for (int i = 0; i < M; ++i) push(A[i]);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << sum << " ";
        if (i+M < N) push(A[i+M]), pop(A[i]);
    }
    cout << endl;
}

// AOJ 3333 - Range Xor Query
// Binary Trie + Sqrt Decomposition
void AOJ3333() {
    const int width = 2000;
    int N, Q;
    cin >> N >> Q;
    vector<int> A(N);
    vector<BinaryTrie<int,12>> vbt(N/width + 1);
    for (int i = 0; i < N; ++i) {
        cin >> A[i];
        vbt[i / width].insert(A[i]);
    }
    
    while (Q--) {
        int type;
        cin >> type;
        if (type == 1) {
            int k, x;
            cin >> k >> x;
            --k;
            vbt[k / width].erase(A[k]);
            A[k] ^= x;
            vbt[k / width].insert(A[k]);
        } else {
            int l, r, x;
            cin >> l >> r >> x;
            --l;
            int res = 1<<29;
            int iter = l;
            while (iter < r && iter % width != 0) {
                res = min(res, A[iter] ^ x);
                ++iter;
            }
            while (iter + width <= r) {
                // ^x をしたときの最小値 (最後に ^x が必要)
                int tmp = vbt[iter / width].get_min(x) ^ x;
                res = min(res, tmp);
                iter += width;
            }
            while (iter < r) {
                res = min(res, A[iter] ^ x);
                ++iter;
            }
            cout << res << endl;
        }
    }
}

// ABC 451 G - Minimum XOR Walk (for insert(), count_lower(add))
// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};
// graph class
template<class T = long long> struct Graph {
    int V, E;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, int m = 0) : V(n), E(m), list(n), reversed_list(n), id(n) { }
    void init(int n = 0, int m = 0) {
        V = n, E = m;
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
        id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;

    // getters
    vector<Edge<T>> &operator [] (int i) { return list[i]; }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
    const void clear() { V = 0; list.clear(); }
    const void resize(int n) { V = n; list.resize(n); }
    Edge<T> &get_edge(int u, int v) {
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<T> &get_edge(int u, int v) const {
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        id[from][to] = (int)list[from].size(), list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    void add_bidirected_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        id[from][to] = (int)list[from].size(), list[from].push_back(Edge(from, to, val));
        reversed_list[from].push_back(Edge(from, to, val));
        if (from != to) {
            id[to][from] = (int)list[to].size(), list[to].push_back(Edge(to, from, val));
            reversed_list[to].push_back(Edge(to, from, val));
        }
    }

    // input / output
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.E; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }
    friend ostream &operator << (ostream &os, const Graph &G) {
        os << endl;
        for (int i = 0; i < G.size(); ++i) {
            os << i << " -> ";
            for (int j = 0; j < G[i].size(); j++) {
                if (j) os << ", ";
                os << G[i][j].to << "(" << G[i][j].val << ")";
            }
            os << endl;
        }
        return os;
    }
};
// cycle basis
// depth[v] := min value of root-v path
// return: cycle basis
template<class T> vector<T> cycle_basis(const Graph<T> &G, vector<T> &depth, int root = 0) {
    int N = (int)G.size();
    depth.assign(N, -1);
    auto rec = [&](auto &&rec, int v) -> void {
        for (auto e : G[v]) {
            if (depth[e.to] != -1) continue;
            depth[e.to] = depth[v] ^ e.val;
            rec(rec, e.to);
        }
    };
    depth[root] = 0;
    rec(rec, root);
    vector<T> base;
    for (int v = 0; v < N; v++) {
        for (auto e : G[v]) {
            T val = depth[e.from] ^ depth[e.to] ^ e.val;
            for (auto b : base) val = min(val, val ^ b);
            if (val) base.emplace_back(val);
        }
    }
    for (auto &val : depth) for (auto b : base) val = min(val, val ^ b);
    return base;
}
void ABC_451_G() {
    int T;
    cin >> T;
    while (T--) {
        long long N, M, K, res = 0;
        cin >> N >> M >> K;
        Graph<long long> G(N);
        for (int i = 0; i < M; i++) {
            long long U, V, W;
            cin >> U >> V >> W, U--, V--;
            G.add_bidirected_edge(U, V, W);
        }
        vector<long long> depth(N, -1);
        auto base = cycle_basis(G, depth, 0);
        BinaryTrie<long long, 50> bt;
        for (int v = 0; v < N; v++) {
            res += bt.count_lower(K+1, depth[v]);
            bt.insert(depth[v]);
        }
        cout << res << endl;
    }
}


int main() {
    YosupoSetXorMin();
    //ARC033_C();
    //AtCoderDwango5Honsen_B();
    //ABC281_E();
    //AOJ3333();
    //ABC_451_G();
}
