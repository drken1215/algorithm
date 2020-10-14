//
// 円と円の共通部分の面積
//
// verified:
//
//   AOJ Course CGL_7_I: Area of Intersection between Two Circles
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_7_I&lang=ja
//
//   Educational Codeforces Round 2 D. Area of Two Circles' Intersection
//     https://codeforces.com/contest/600/problem/D
//


#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


////////////////////////////
// 基本要素 (点, 線分, 円)
////////////////////////////

using DD = long double;
const DD INF = 1LL<<60;      // to be set appropriately
const DD EPS = 1e-6;        // to be set appropriately
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


// 円と円の共通部分の面積
DD calc(const Circle &p, const Circle &q) {
    DD d = abs(p - q);
    if (d >= p.r + q.r - EPS) return 0;
    else if (d <= abs(p.r - q.r) + EPS) return min(p.r, q.r) * min(p.r, q.r) * PI;
    DD pcos = (p.r*p.r + d*d - q.r*q.r) / (p.r*d*2);
    DD pang = acosl(pcos);
    DD parea = p.r*p.r*pang - p.r*p.r*sin(pang*2)/2;
    DD qcos = (q.r*q.r + d*d - p.r*p.r) / (q.r*d*2);
    DD qang = acosl(qcos);
    DD qarea = q.r*q.r*qang - q.r*q.r*sin(qang*2)/2;
    return parea + qarea;
}



///////////////////////////////////////////////
// solver
///////////////////////////////////////////////

int main() {
    Circle p, q;
    while (cin >> p.x >> p.y >> p.r >> q.x >> q.y >> q.r) {
        cout << fixed << setprecision(15) << calc(p, q) << endl;
    }
}
