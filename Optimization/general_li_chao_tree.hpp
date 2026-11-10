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
    
    // inner values
    vector<unique_ptr<Node[]>> blocks;  // chunk-based memory pool
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

// General Li-Chao Tree (Static ver.) for Pairwise Single Crossing Functions
template<class VAL, class FUNC, class X = long long> struct GeneralLiChaoTree {
    // inner values
    int N;
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<X> xs;
    vector<int> ids;
    vector<FUNC> funcs;

    // constructors
    explicit GeneralLiChaoTree(vector<X> xs_, VAL inf = numeric_limits<VAL>::max() / 2, bool sorted = true) 
    : xs(std::move(xs_)), INF(inf) {
        if (!sorted) {
            sort(xs.begin(), xs.end());
            xs.erase(unique(xs.begin(), xs.end()), xs.end());
        }
        N = (int)xs.size();
        ids.assign(max(1, N * 4), -1);
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
        if (N) insert_at(1, 0, N, id);
        return id;
    }

    // insert_segment(f, l, r), in O(log^2 N)
    int insert_segment(FUNC f, X l, X r) {
        assert(l <= r);
        int id = (int)funcs.size();
        funcs.push_back(std::move(f));
        if (!N || l == r) return id;
        int L = int(lower_bound(xs.begin(), xs.end(), l) - xs.begin());
        int R = int(lower_bound(xs.begin(), xs.end(), r) - xs.begin());
        if (L < R) insert_segment_at(1, 0, N, L, R, id);
        return id;
    }

    // query(x), x must be in xs, in O(log N)
    pair<VAL, int> query(X x) const {
        auto it = lower_bound(xs.begin(), xs.end(), x);
        assert(it != xs.end() && *it == x);
        int pos = int(it - xs.begin());
        VAL best = INF;
        int best_id = -1;
        int v = 1, left = 0, right = N;
        while (true) {
            int fid = ids[v];
            if (fid != -1) {
                VAL y = eval(fid, x);
                if (best_id == -1 || y < best || (y == best && fid < best_id)) {
                    best = y, best_id = fid;
                }
            }
            if (right - left == 1) break;
            int mid = (left + right) / 2;
            if (pos < mid) v *= 2, right = mid;
            else v = v * 2 + 1, left = mid;
        }
        return {best, best_id};
    }

private:
    // inner insert
    void insert_at(int v, int left, int right, int fid) {
        while (true) {
            int gid = ids[v];
            if (gid == -1) {
                ids[v] = fid;
                return;
            }
            bool left_better = better(fid, gid, xs[left]);
            bool right_better = better(fid, gid, xs[right - 1]);
            if (!left_better && !right_better) return;
            if (left_better && right_better) {
                ids[v] = fid;
                return;
            }
            int mid = (left + right) / 2;
            bool mid_better = better(fid, gid, xs[mid]);
            if (mid_better) ids[v] = fid, fid = gid;
            if (left_better != mid_better) v *= 2, right = mid;
            else v = v * 2 + 1, left = mid;
        }
    }

    // inner insert_segment
    void insert_segment_at(int v, int left, int right, int l, int r, int fid) {
        if (r <= left || right <= l) return;
        if (l <= left && right <= r) {
            insert_at(v, left, right, fid);
            return;
        }
        int mid = left + (right - left) / 2;
        if (l < mid) insert_segment_at(v * 2, left, mid, l, r, fid);
        if (mid < r) insert_segment_at(v * 2 + 1, mid, right, l, r, fid);
    }
};
