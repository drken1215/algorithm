#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// General Li-Chao Tree (Dynamic ver.) for Pairwise Single Crossing Functions
template<class VAL, class FUNC, class X = long long> struct DynamicGeneralLiChaoTree {
    struct Node {
        int id = -1;  // function id
        Node *left = nullptr, *right = nullptr;
    };
    
    // inner values using chunk-based memory pool
    vector<unique_ptr<Node[]>> blocks;
    vector<FUNC> funcs;
    Node *root = nullptr;
    X XL, XR;
    VAL INF = numeric_limits<VAL>::max() / 2;
    
    // constructors
    explicit DynamicGeneralLiChaoTree(X xl, X xr, VAL inf = numeric_limits<VAL>::max() / 2)
    : XL(xl), XR(xr), INF(inf) {
        assert(xl < xr);
    }

    // basic functions
    void reserve_functions(size_t n) {
        funcs.reserve(n);
    }
    VAL eval(int id, X x) const {
        return funcs[id](x);
    }
    bool better(int fid, int gid, X x) const {
        VAL a = eval(fid, x), b = eval(gid, x);
        return (a < b) || (a == b && fid < gid);
    }
    void clear() {
        funcs.clear();
        blocks.clear();
        used = BLOCK;
        root = nullptr;
    }

    // insert(f), in O(log X)
    int insert(FUNC f) {
        int id = (int)funcs.size();
        funcs.push_back(std::move(f));
        insert_at(root, XL, XR, id);
        return id;
    }

    // insert_segment(l, r, f), in O(log^2 X)
    int insert_segment(FUNC f, X l, X r) {
        assert(XL <= l && l <= r && r <= XR);
        int id = (int)funcs.size();
        funcs.push_back(std::move(f));
        if (l != r) insert_segment_at(root, XL, XR, l, r, id);
        return id;
    }

    // query(x), in O(log X)
    pair<VAL, int> query(X x) const {
        assert(XL <= x && x < XR);
        VAL best = INF;
        int best_id = -1;
        Node *p = root;
        X left = XL, right = XR;
        while (p) {
            if (p->id != -1) {
                int fid = p->id;
                VAL y = eval(fid, x);
                if (best_id == -1 || y < best || (y == best && fid < best_id)) {
                    best = y, best_id = fid;
                }
            }
            if (right - left == 1) break;
            X mid = (left + right) / 2;
            if (x < mid) p = p->left, right = mid;
            else p = p->right, left = mid;
        }
        return {best, best_id};
    }

private:
    // chunk-based node allocation
    static constexpr int BLOCK = 4096;
    int used = BLOCK;
    Node *new_node() {
        if (used == BLOCK) {
            blocks.emplace_back(make_unique<Node[]>(BLOCK));
            used = 0;
        }
        return &blocks.back()[used++];
    }

    // inner insert
    void insert_at(Node *&p, X left, X right, int fid) {
        Node **cur = &p;
        while (true) {
            if (!*cur) *cur = new_node();
            Node *v = *cur;
            int gid = v->id;
            if (gid == -1) {
                v->id = fid;
                return;
            }
            bool left_better = better(fid, gid, left), right_better = better(fid, gid, right - 1);
            if (!left_better && !right_better) return;
            if (left_better && right_better) {
                v->id = fid;
                return;
            }
            X mid = (left + right) / 2;
            bool mid_better = better(fid, gid, mid);
            if (mid_better) v->id = fid, fid = gid;
            if (left_better != mid_better) cur = &v->left, right = mid;
            else cur = &v->right, left = mid;
        }
    }

    // inner insert_segment
    void insert_segment_at(Node *&p, X left, X right, X l, X r, int fid) {
        if (r <= left || right <= l) return;
        if (l <= left && right <= r) {
            insert_at(p, left, right, fid);
            return;
        }
        if (!p) p = new_node();
        X mid = (left + right) / 2;
        if (l < mid) insert_segment_at(p->left, left, mid, l, r, fid);
        if (mid < r) insert_segment_at(p->right, mid, right, l, r, fid);
    }
};

// General Li-Chao Tree (Coordinate Compression ver.) for Pairwise Single Crossing Functions
template<class VAL, class FUNC, class X = long long> struct GeneralLiChaoTree {
    // inner values
    int N, XL, SIZ;
    bool compressed = true;
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<X> xs;
    vector<int> ids;
    vector<FUNC> funcs;

    // constructors
    explicit GeneralLiChaoTree(vector<X> xs_, VAL inf = numeric_limits<VAL>::max() / 2, bool sorted = true) 
    : N((int)xs_.size()), XL(0), SIZ(1), compressed(true), INF(inf), xs(std::move(xs_)) {
        if (!sorted) {
            sort(xs.begin(), xs.end());
            xs.erase(unique(xs.begin(), xs.end()), xs.end());
        }
        while (SIZ < N) SIZ <<= 1;
        ids.assign(SIZ << 1, -1);
    }
    explicit GeneralLiChaoTree(int xl, int xr, VAL inf = numeric_limits<VAL>::max() / 2) 
    : N(xr - xl), XL(xl), SIZ(1), compressed(false), INF(inf) {
        assert(N >= 0);
        while (SIZ < N) SIZ <<= 1;
        ids.assign(SIZ << 1, -1);
    }

    // basic functions
    void reserve_functions(size_t n) {
        funcs.reserve(n);
    }
    VAL eval(int id, X x) const {
        return funcs[id](x);
    }
    bool better(int fid, int gid, X x) const {
        VAL a = eval(fid, x), b = eval(gid, x);
        return (a < b) || (a == b && fid < gid);
    }
    void clear() {
        funcs.clear();
        fill(ids.begin(), ids.end(), -1);
    }

    // insert(f), in O(log N)
    int insert(FUNC f) {
        int id = (int)funcs.size();
        funcs.push_back(std::move(f));
        if (N) insert_at(1, 0, SIZ, id);
        return id;
    }

    // insert_segment(f, l, r), in O(log^2 N)
    int insert_segment(FUNC f, X l, X r) {
        assert(l <= r);
        int id = (int)funcs.size();
        funcs.push_back(std::move(f));
        if (N == 0 || l == r) return id;
        int L = lower_bound_index(l);
        int R = lower_bound_index(r);
        if (L == R) return id;
        int a = L + SIZ, b = R + SIZ;
        while (a < b) {
            if (a & 1) {
                int v = a++;
                int depth = 31 - __builtin_clz((unsigned)v);
                int length = SIZ >> depth;
                int left = (v - (1 << depth)) * length;
                insert_at(v, left, left + length, id);
            }
            if (b & 1) {
                int v = --b;
                int depth = 31 - __builtin_clz((unsigned)v);
                int length = SIZ >> depth;
                int left = (v - (1 << depth)) * length;
                insert_at(v, left, left + length, id);
            }
            a >>= 1, b >>= 1;
        }
        return id;
    }

    // query(x), x must be in xs, in O(log N)
    pair<VAL, int> query(X x) const {
        int pos = 0;
        if (compressed) {
            auto it = lower_bound(xs.begin(), xs.end(), x);
            assert(it != xs.end() && *it == x);
            pos = int(it - xs.begin());
        } else {
            assert(X(XL) <= x && x < X(XL + N));
            pos = int(x - XL);
        }
        VAL best = INF;
        int best_id = -1;
        int v = pos + SIZ;
        while (v) {
            int fid = ids[v];
            if (fid != -1) {
                VAL y = eval(fid, x);
                if (best_id == -1 || y < best || (y == best && fid < best_id)) {
                    best = y, best_id = fid;
                }
            }
            v >>= 1;
        }
        return {best, best_id};
    }

private:
    // coordinate compression
    int lower_bound_index(X x) const {
        if (compressed) {
            return int(std::lower_bound(xs.begin(), xs.end(), x) - xs.begin());
        } else {
            if (x <= X(XL)) return 0;
            else if (x >= X(XL) + X(N)) return N;
            else return int(x - X(XL));
        }
    }

    // inner insert
    void insert_at(int v, int left, int right, int fid) {
        auto coord = [&](int i) -> X { 
            i = min(i, N - 1);
            return (compressed ? xs[i] : X(i) + X(XL));
        };
        while (true) {
            int gid = ids[v];
            if (gid == -1) {
                ids[v] = fid;
                return;
            }
            bool left_better = better(fid, gid, coord(left));
            bool right_better = better(fid, gid, coord(right - 1));
            if (!left_better && !right_better) return;
            if (left_better && right_better) {
                ids[v] = fid;
                return;
            }
            X mid = (left + right) / 2;
            bool mid_better = better(fid, gid, coord(mid));
            if (mid_better) ids[v] = fid, fid = gid;
            if (left_better != mid_better) {
                v <<= 1;
                right = mid;
            } else {
                v = v << 1 | 1;
                left = mid;
            }
        }
    }
};
