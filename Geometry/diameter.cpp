//
// 凸多角形の直径, O(n)
//
// verified:
//   AOJ Course CGL_4_B Diameter of a Convex Polygon
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_B&lang=ja
//


//------------------------------//
// 基本要素 (点, 線分, 円)
//------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;


using DD = double;
const DD INF = 1LL<<60;      // to be set appropriately
const DD EPS = 1e-10;        // to be set appropriately
const DD PI = acosl(-1.0);
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

/* Circle */
struct Circle : Point {
    DD r;
    Circle(Point p = Point(0.0, 0.0), DD r = 0.0) : Point(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle &c) {return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';}
};


//------------------------------//
// 多角形
//------------------------------//

// 一直線上の3点を含めない
vector<Point> ConvexHull(vector<Point> &ps) {
    int n = (int)ps.size();
    vector<Point> res(2*n);
    sort(ps.begin(), ps.end());
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (k >= 2) {
            while (cross(res[k-1] - res[k-2], ps[i] - res[k-2]) < EPS) {
                --k;
                if (k < 2) break;
            }
        }
        res[k] = ps[i]; ++k;
    }
    int t = k+1;
    for (int i = n-2; i >= 0; --i) {
        if (k >= t) {
            while (cross(res[k-1] - res[k-2], ps[i] - res[k-2]) < EPS) {
                --k;
                if (k < t) break;
            }
        }
        res[k] = ps[i]; ++k;
    }
    res.resize(k-1);
    return res;
}

// 凸多角形 pol の直径を求める
DD CalcDiameter(const vector<Point> &pol) {
    int n = (int)pol.size();
    int si = 0, sj = 0;
    for (int i = 0; i < n; ++i) {
        if (pol[i].y < pol[si].y) si = i;
        if (pol[i].y > pol[sj].y) sj = i;
    }
    
    DD res = 0;
    int i = si, j = sj;
    do {
        res = max(res, abs(pol[i] - pol[j]));
        if (cross(pol[(i+1)%n] - pol[i], pol[(j+1)%n] - pol[j]) >= 0)
            j = (j+1)%n;
        else
            i = (i+1)%n;
    } while (i != si || j != sj);
    
    return res;
}


//------------------------------//
// Examples
//------------------------------//

void AOJ_CGL_4_B() {
    int n;
    cin >> n;
    vector<Point> pol(n);
    for (int i = 0; i < n; ++i) cin >> pol[i].x >> pol[i].y;
    cout << fixed << setprecision(10) << CalcDiameter(pol) << endl;
}

int main() {
    AOJ_CGL_4_B();
}
