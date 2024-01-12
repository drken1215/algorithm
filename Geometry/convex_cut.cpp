//
// 凸多角形の切断, O(n)
//
// verified:
//   AOJ Course CGL_4_C Convex Polygon - Convex Cut
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_C&lang=jp
//

/*
    問題例
    ・AOJ 2385
    ・AOJ 2160
*/


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

// 多角形の面積
DD Area(const vector<Point> &pol) {
    DD res = 0.0;
    for (int i = 0; i < pol.size(); ++i) {
        res += cross(pol[i], pol[(i+1)%pol.size()]);
    }
    return res/2.0L;
}

// convex cut
int ccw_for_convexcut(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
vector<Point> crosspoint_for_convexcut(const Line &l, const Line &m) {
    vector<Point> res;
    DD d = cross(m[1] - m[0], l[1] - l[0]);
    if (abs(d) < EPS) return vector<Point>();
    res.push_back(l[0] + (l[1] - l[0]) * cross(m[1] - m[0], m[1] - l[0]) / d);
    return res;
}
vector<Point> ConvexCut(const vector<Point> &pol, const Line &l) {
    vector<Point> res;
    for (int i = 0; i < pol.size(); ++i) {
        Point p = pol[i], q = pol[(i+1)%pol.size()];
        if (ccw_for_convexcut(l[0], l[1], p) != -1) {
            if (res.size() == 0) res.push_back(p);
            else if (!eq(p, res[res.size()-1])) res.push_back(p);
        }
        if (ccw_for_convexcut(l[0], l[1], p) * ccw_for_convexcut(l[0], l[1], q) < 0) {
            vector<Point> temp = crosspoint_for_convexcut(Line(p, q), l);
            if (temp.size() == 0) continue;
            else if (res.size() == 0) res.push_back(temp[0]);
            else if (!eq(temp[0], res[res.size()-1])) res.push_back(temp[0]);
        }
    }
    return res;
}

// Voronoi-diagram (今回は使わない)
Line bisector(const Point &p, const Point &q) {
    Point c = (p + q) / 2.0L;
    Point v = (q - p) * Point(0.0L, 1.0L);
    v = v / abs(v);
    return Line(c - v, c + v);
}
vector<Point> Voronoi(vector<Point> pol, const vector<Point> &ps, int ind) {
    for (int i = 0; i < ps.size(); ++i) {
        if (i == ind) continue;
        Line l = bisector(ps[ind], ps[i]);
        pol = ConvexCut(pol, l);
    }
    return pol;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int n; cin >> n;
    vector<Point> pol(n);
    for (int i = 0; i < n; ++i) cin >> pol[i].x >> pol[i].y;
    int Q; cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        Point x, y; cin >> x.x >> x.y >> y.x >> y.y;
        Line l(x, y);
        auto cutted = ConvexCut(pol, l);
        cout << fixed << setprecision(10) << Area(cutted) << endl;
    }
}
