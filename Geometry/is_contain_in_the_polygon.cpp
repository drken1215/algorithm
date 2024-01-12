//
// 点が多角形 (凸とは限らない) に含まれるかどうか, O(n)
//
// verified:
//   AOJ Course CGL_3_C Polygon - Polygon-Point Containment
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_3_C&lang=jp
//


#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;


//------------------------------//
// 基本要素 (点, 線分, 円)
//------------------------------//

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
// 点と多角形の包含
//------------------------------//

// 2: in, 1: on, 0: out
int is_contain(const vector<Point> &pol, const Point &p) {
    int n = (int)pol.size();
    int isin = 0;
    for (int i = 0; i < n; ++i) {
        Point a = pol[i] - p, b = pol[(i+1)%n] - p;
        if (a.y > b.y) swap(a, b);
        if (a.y <= 0 && b.y > 0) if (cross(a, b) < 0) isin = 1-isin;
        if (cross(a, b) == 0 && dot(a, b) <= 0) return 1;
    }
    if (isin) return 2;
    else return 0;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int n;
    cin >> n;
    vector<Point> pol(n);
    for (int i = 0; i < n; ++i) cin >> pol[i].x >> pol[i].y;
    int Q; cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        Point p; cin >> p.x >> p.y;
        cout << is_contain(pol, p) << endl;
    }
}
