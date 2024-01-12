//
// 最近点対
//
// verified:
//   AOJ Course CGL_5_A Point Set - Closest Pair
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_5_A&lang=jp
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


// 最近点対
bool compare_y(Point a, Point b) { return a.y < b.y; }
DD DivideAndConqur(vector<Point>::iterator it, int n) {
    if (n <= 1) return INF;
    int m = n/2;
    DD x = it[m].x;
    DD d = min(DivideAndConqur(it, m), DivideAndConqur(it+m, n-m));
    inplace_merge(it, it+m, it+n, compare_y);
    
    vector<Point> vec;
    for (int i = 0; i < n; ++i) {
        if (fabs(it[i].x - x) >= d) continue;
        for (int j = 0; j < vec.size(); ++j) {
            DD dx = it[i].x - vec[vec.size()-j-1].x;
            DD dy = it[i].y - vec[vec.size()-j-1].y;
            if (dy >= d) break;
            d = min(d, sqrt(dx*dx+dy*dy));
        }
        vec.push_back(it[i]);
    }
    return d;
}
DD Closet(vector<Point> ps) {
    int n = (int)ps.size();
    sort(ps.begin(), ps.end());
    return DivideAndConqur(ps.begin(), n);
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int N; cin >> N;
    vector<Point> ps(N);
    for (int i = 0; i < N; ++i) cin >> ps[i].x >> ps[i].y;
    cout << fixed << setprecision(10) << Closet(ps) << endl;
}
