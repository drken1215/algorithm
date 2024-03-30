//
// 二次元セグメント木
//
// verified:
//   ABC 347 F - Non-overlapping Squares
//     https://atcoder.jp/contests/abc347/tasks/abc347_f
//


#include <bits/stdc++.h>
using namespace std;


// Segment Tree 2D
template<class Monoid> struct SegmentTree2D {
    using Func = function<Monoid(Monoid, Monoid)>;
    
    // core member
    int H, W;
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int log_H, offset_H, log_W, offset_W;
    vector<Monoid> dat;
    
    // constructor
    SegmentTree2D() {}
    SegmentTree2D(int h, int w, const Func &op, const Monoid &identity) {
        init(h, w, op, identity);
    }
    SegmentTree2D(const vector<vector<Monoid>> &v, const Func &op, const Monoid &identity) {
        init(v, op, identity);
    }
    void init(int h, int w, const Func &op, const Monoid &identity) {
        H = h, W = w;
        OP = op;
        IDENTITY = identity;
        log_H = 0, offset_H = 1, log_W = 0, offset_W = 1;
        while (offset_H < H) ++log_H, offset_H <<= 1;
        while (offset_W < W) ++log_W, offset_W <<= 1;
        dat.assign(offset_H * offset_W * 4, IDENTITY);
    }
    void init(const vector<vector<Monoid>> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), (int)v[0].size(), op, identity);
        build(v);
    }
    int id(int x, int y) const {
        return x * 2 * offset_W + y;
    }
    void build(const vector<vector<Monoid>> &v) {
        assert(H == (int)v.size()), assert(W == (int)v[0].size());
        for (int x = 0; x < H; ++x) {
            for (int y = 0; y < W; ++y) {
                dat[id(x + offset_H, y + offset_W)] = v[x][y];
            }
        }
        for (int y = offset_W; y < offset_W * 2; ++y) {
            for (int x = offset_H - 1; x; --x) {
                dat[id(x, y)] = OP(dat[id(x * 2, y)], dat[id(x * 2 + 1, y)]);
            }
        }
        for (int x = 0; x < offset_H * 2; ++x) {
            for (int y = offset_W - 1; y; --y) {
                dat[id(x, y)] = OP(dat[id(x, y * 2)], dat[id(x, y * 2 + 1)]);
            }
        }
    }
    Monoid operator () (int x, int y) const {
        return dat[id(x + offset_H, y + offset_W)];
    }
    
    // update
    void set(int x, int y, const Monoid &v) {
        assert(x >= 0 && x < H && y >= 0 && y < W);
        x += offset_H, y += offset_W;
        dat[id(x, y)] = v;
        for (int i = x >> 1; i; i >>= 1) {
            dat[id(i, y)] = OP(dat[id(i * 2, y)], dat[id(i * 2 + 1, y)]);
        }
        for (; x; x >>= 1) {
            for (int j = y >> 1; j; j >>= 1) {
                dat[id(x, j)] = OP(dat[id(x, j * 2)], dat[id(x, j * 2 + 1)]);
            }
        }
    }
    
    // prod
    Monoid inner_prod(int x, int yl, int yr) {
        Monoid res = IDENTITY;
        for (; yl < yr; yl >>= 1, yr >>= 1) {
            if (yl & 1) res = OP(res, dat[id(x, yl++)]);
            if (yr & 1) res = OP(dat[id(x, --yr)], res);
        }
        return res;
    }
    Monoid prod(int xl, int xr, int yl, int yr) {
        assert(0 <= xl && xl <= xr && xr <= H);
        assert(0 <= yl && yl <= yr && yr <= W);
        Monoid res = IDENTITY;
        xl += offset_H, xr += offset_H, yl += offset_W, yr += offset_W;
        for (; xl < xr; xl >>= 1, xr >>= 1) {
            if (xl & 1) res = OP(res, inner_prod(xl++, yl, yr));
            if (xr & 1) res = OP(inner_prod(--xr, yl, yr), res);
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream &s, const SegmentTree2D &seg) {
        for (int x = 0; x < seg.H; ++x) {
            for (int y = 0; y < seg.W; ++y) {
                s << seg(x, y) << " ";
            }
            s << endl;
        }
        return s;
    }
};



//------------------------------//
// Examples
//------------------------------//

// ABC 347 F - Non-overlapping Squares
void ABC_347_F() {
    const long long INF = 1LL<<55;
    long long res = 0;
    int N, M;
    cin >> N >> M;
    vector A(N, vector(N, 0LL)), S(N+1, vector(N+1, 0LL)), v(N-M+1, vector(N-M+1, -INF));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cin >> A[i][j];
            S[i+1][j+1] = S[i+1][j] + S[i][j+1] - S[i][j] + A[i][j];
        }
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i+M <= N && j+M <= N) {
                v[i][j] = S[i+M][j+M] - S[i+M][j] - S[i][j+M] + S[i][j];
            }
        }
    }
    auto v2 = v;
    int HW = (int)v.size();
    
    auto op = [](long long a, long long b) {return max(a, b);};
    {
        SegmentTree2D<long long> seg(v, op, -INF);
        for (int i = 0; i < HW; ++i) {
            for (int j = i+M; j+M < HW; ++j) {
                res = max(res, seg.prod(0, HW, i, i+1) + seg.prod(0, HW, j, j+1)
                          + seg.prod(0, HW, j+M, HW));
            }
        }
        for (int i = 0; i < HW; ++i) {
            for (int j = i+M; j+M < HW; ++j) {
                res = max(res, seg.prod(i, i+1, 0, HW) + seg.prod(j, j+1, 0, HW)
                          + seg.prod(j+M, HW, 0, HW));
            }
        }
    }

    auto calc = [&]() {
        SegmentTree2D<long long> seg(v, op, -INF);
        long long top = 0;
        for (int x = 1; x <= HW; ++x) {
            for (int y = 0; y+M < HW; ++y) {
                top = max(top, seg.prod(0, x, y, y+1) + seg.prod(0, x, y+M, HW));
            }
            if (x+M-1 < HW) res = max(res, top + seg.prod(x+M-1, HW, 0, HW));
        }
    };
    
    calc();
    reverse(v.begin(), v.end());
    calc();
    for (int i = 0; i < HW; ++i) for (int j = 0; j < HW; ++j) v[i][j] = v2[j][i];
    calc();
    reverse(v.begin(), v.end());
    calc();

    cout << res << endl;
}

int main() {
    ABC_347_F();
}

