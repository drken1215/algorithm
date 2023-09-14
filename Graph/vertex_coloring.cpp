//
// 彩色数を求める O(N2^N) のアルゴリズム
//
// cf.
//   https://drken1215.hatenablog.com/entry/2019/01/16/030000
//
// verified
//   AOJ 2136 Webby Subway
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2136
//


#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


// 彩色数
long long modpow(long long a, long long n, long long MOD) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % MOD;
        a = a * a % MOD;
        n >>= 1;
    }
    return res;
}
int chromatic_number(const vector<vector<int> > &G) {
    const int MOD = 1000000007;
    int n = (int)G.size();
    vector<int> neighbor(n, 0);
    for (int i = 0; i < n; ++i) {
        int S = (1<<i);
        for (int j = 0; j < n; ++j)
            if (G[i][j])
                S |= (1<<j);
        neighbor[i] = S;
    }
    
    // I[S] := #. of inndepndet subset of S
    vector<int> I(1<<n);
    I[0] = 1;
    for (int S = 1; S < (1<<n); ++S) {
        int v = __builtin_ctz(S);
        I[S] = I[S & ~(1<<v)] + I[S & ~neighbor[v]];
    }
    int low = 0, high = n;
    while (high - low > 1) {
        int mid = (low + high) >> 1;
        
        // g[S] := #. of "k independent sets" which cover S just
        // f[S] := #. of "k independent sets" which consist of subseta of S
        // then
        //   f[S] = sum_{T in S} g(T)
        //   g[S] = sum_{T in S} (-1)^(|S|-|T|)f[T]
        long long g = 0;
        for (int S = 0; S < (1<<n); ++S) {
            if ((n - __builtin_popcount(S)) & 1) g -= modpow(I[S], mid, MOD);
            else g += modpow(I[S], mid, MOD);
            g = (g % MOD + MOD) % MOD;
        }
        if (g != 0) high = mid;
        else low = mid;
    }
    return high;
}



////////////////////////////
// solver
////////////////////////////

using DD = double;
const DD INF = 1LL<<60;      // to be set appropriately
const DD EPS = 1e-10;        // to be set appropriately
const DD PI = acos(-1.0);
DD torad(int deg) {return (DD)(deg) * PI / 180;}
DD todeg(DD ang) {return ang * 180 / PI;}

/* Point */
struct Point {
    DD x, y;
    Point(DD x = 0.0, DD y = 0.0) : x(x), y(y) {}
    friend ostream& operator << (ostream &s, const Point &p) {return s << '(' << p.x << ", " << p.y << ')';}
};
inline Point operator + (const Point &p, const Point &q) {return Point(p.x + q.x, p.y + q.y);}
inline Point operator - (const Point &p, const Point &q) {return Point(p.x - q.x, p.y - q.y);}
inline Point operator * (const Point &p, DD a) {return Point(p.x * a, p.y * a);}
inline Point operator * (DD a, const Point &p) {return Point(a * p.x, a * p.y);}
inline Point operator * (const Point &p, const Point &q) {return Point(p.x * q.x - p.y * q.y, p.x * q.y + p.y * q.x);}
inline Point operator / (const Point &p, DD a) {return Point(p.x / a, p.y / a);}
inline Point conj(const Point &p) {return Point(p.x, -p.y);}
inline Point rot(const Point &p, DD ang) {return Point(cos(ang) * p.x - sin(ang) * p.y, sin(ang) * p.x + cos(ang) * p.y);}
inline Point rot90(const Point &p) {return Point(-p.y, p.x);}
inline DD cross(const Point &p, const Point &q) {return p.x * q.y - p.y * q.x;}
inline DD dot(const Point &p, const Point &q) {return p.x * q.x + p.y * q.y;}
inline DD norm(const Point &p) {return dot(p, p);}
inline DD abs(const Point &p) {return sqrt(dot(p, p));}
inline DD amp(const Point &p) {DD res = atan2(p.y, p.x); if (res < 0) res += PI*2; return res;}
inline bool eq(const Point &p, const Point &q) {return abs(p - q) < EPS;}
inline bool operator < (const Point &p, const Point &q) {return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);}
inline bool operator > (const Point &p, const Point &q) {return (abs(p.x - q.x) > EPS ? p.x > q.x : p.y > q.y);}
inline Point operator / (const Point &p, const Point &q) {return p * conj(q) / norm(q);}

/* Line */
struct Line : vector<Point> {
    Line(Point a = Point(0.0, 0.0), Point b = Point(0.0, 0.0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line &l) {return s << '{' << l[0] << ", " << l[1] << '}';}
};

int ccw_for_dis(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
bool isinterPS(const Point &p, const Line &s) {
    return (ccw_for_dis(s[0], s[1], p) == 0);
}
bool isinterSS(const Line &s, const Line &t) {
    if (eq(s[0], s[1])) return isinterPS(s[0], t);
    if (eq(t[0], t[1])) return isinterPS(t[0], s);
    return (ccw_for_dis(s[0], s[1], t[0]) * ccw_for_dis(s[0], s[1], t[1]) <= 0 &&
            ccw_for_dis(t[0], t[1], s[0]) * ccw_for_dis(t[0], t[1], s[1]) <= 0);
}

int main() {
    int N;
    while (cin >> N, N) {
        vector<vector<int> > G(N, vector<int>(N, 0));
        vector<vector<Line> > lines(N);
        for (int i = 0; i < N; ++i) {
            int num; cin >> num;
            double x, y; cin >> x >> y;
            for (int j = 1; j < num; ++j) {
                double nx, ny; cin >> nx >> ny;
                Line l(Point(x, y), Point(nx, ny));
                lines[i].push_back(l);
                x = nx, y = ny;
            }
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                bool ok = false;
                for (auto li : lines[i]) {
                    for (auto lj : lines[j]) {
                        if (isinterSS(li, lj)) ok = true;
                    }
                }
                if (ok) G[i][j] = G[j][i] = true;
            }
        }
        cout << chromatic_number(G) << endl;
    }
}
