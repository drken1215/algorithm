#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;


// Dynamic Li Chao Tree for Pairwise Single Crossing Functions
/*
    Func: Pairwise Single Crossing 性を仮定
    注意点：i < j なる領域だけでなく、全領域での Pairwise Single Crossing 性が必要！！！

    x 座標のとりうる値の範囲: [XL, XR)
*/
template<class VAL, class FUNC, class X = long long> struct DynamicGeneralLiChaoTree {
    using Func = function<VAL(X)>;
    struct IndexedFunc {
        Func f;
        int id;
        VAL eval(X x) const { return f(x); };
    };
    struct Node {
        optional<IndexedFunc> func;
        unique_ptr<Node> left, right;
    };
    
    // inner values
    X XL, XR;
    VAL INF = numeric_limits<VAL>::max() / 2;
    int next_id = 0;
    unique_ptr<Node> root;
    
    // constructors
    DynamicGeneralLiChaoTree() {}
    DynamicGeneralLiChaoTree(X xl, X xr, VAL inf = numeric_limits<VAL>::max() / 2) : XL(xl), XR(xr), INF(inf) {
        assert(xl < xr);
    }
    static bool better(const IndexedFunc &a, const IndexedFunc &b, X x) {
        return make_pair(a.eval(x), a.id) < make_pair(b.eval(x), b.id);
    }
    void clear() {
        root.reset();
        next_id = 0;
    }

    // insert(f), in O(log X)
    void insert(IndexedFunc f, unique_ptr<Node> &p, X left, X right) {
        if (!p) p = make_unique<Node>();
        if (!p->func) {
            p->func.emplace(f);
            return;
        }
        X mid = (left + right) / 2;
        bool left_better = better(f, *p->func, left);
        bool right_better = better(f, *p->func, right - 1);

        if (!left_better && !right_better) return;
        if (left_better && right_better) {
            p->func.emplace(f);
            return;
        }
        bool mid_better = better(f, *p->func, mid);
        if (mid_better) {
            IndexedFunc old = std::move(*p->func);
            p->func.emplace(f);
            if (left_better != mid_better) insert(old, p->left, left, mid);
            else insert(old, p->right, mid, right);
        } else {
            if (left_better != mid_better) insert(f, p->left, left, mid);
            else insert(f, p->right, mid, right);
        }
    }
    int insert(Func f) {
        int id = next_id++;
        insert(IndexedFunc{std::move(f), id}, root, XL, XR);
        return id;
    }

    // insert_segment(l, r, f), in O(log^2 X)
    void insert_segment(const IndexedFunc &f, unique_ptr<Node> &p, X left, X right, X l, X r) {
        if (r <= left || right <= l) return;
        if (l <= left && right <= r) {
            insert(f, p, left, right);
            return;
        }
        if (!p) p = make_unique<Node>();
        X mid = (left + right) / 2;
        insert_segment(f, p->left, left, mid, l, r);
        insert_segment(f, p->right, mid, right, l, r);
    }
    int insert_segment(Func f, X l, X r) {
        assert(XL <= l && l <= r && r <= XR);
        int id = next_id++;
        if (l != r) insert_segment(IndexedFunc{std::move(f), id}, root, XL, XR, l, r);
        return id;
    }

    // query(x), in O(log X)
    pair<VAL, int> query(X x) const {
        assert(XL <= x && x < XR);
        pair<VAL, int> res{INF, -1};
        const Node *p = root.get();
        X left = XL, right = XR;
        while (p) {
            if (p->func) {
                auto cur = make_pair(p->func->eval(x), p->func->id);
                if (res.second == -1 || cur < res) res = cur;
            }
            if (right - left == 1) break;
            X mid = (left + right) / 2;
            if (x < mid) {
                p = p->left.get();
                right = mid;
            } else {
                p = p->right.get();
                left = mid;
            }
        }
        return res;
    }
};

