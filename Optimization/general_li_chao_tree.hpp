#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;


// Dynamic General Li-Chao Tree for Pairwise Single Crossing Functions
template<class VAL, class FUNC, class X = long long> struct DynamicGeneralLiChaoTree {
    struct Node {
        int id = -1;  // function id
        Node *left = nullptr, *right = nullptr;
    };
    
    // inner values using chunk-based memory pool
    static constexpr int BLOCK = 4096;
    vector<unique_ptr<Node[]>> blocks;
    int used = BLOCK;
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
    Node *new_node() {
        if (used == BLOCK) {
            blocks.emplace_back(make_unique<Node[]>(BLOCK));
            used = 0;
        }
        return &blocks.back()[used++];
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
                if (best_id == -1 || y < best || (y == best && fid < best_id)) best = y, best_id = fid;
            }
            if (right - left == 1) break;
            X mid = (left + right) / 2;
            if (x < mid) p = p->left, right = mid;
            else p = p->right, left = mid;
        }
        return {best, best_id};
    }

private:
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
