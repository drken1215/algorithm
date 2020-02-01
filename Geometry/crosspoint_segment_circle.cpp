//
// 円と線分の交点
//
// verified:
//   AOJ 1183 鎖中経路
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1183&lang=jp
//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;


////////////////////////////
// 基本要素 (点, 線分, 円)
////////////////////////////

using DD = double;
const DD INF = 1LL << 60;      // to be set appropriately
const DD EPS = 1e-5;        // to be set appropriately
const DD PI = acosl(-1.0);
DD torad(int deg) { return (DD)(deg)* PI / 180; }
DD todeg(DD ang) { return ang * 180 / PI; }

/* Point */
struct Point {
    DD x, y;
    Point(DD x = 0.0, DD y = 0.0) : x(x), y(y) {}
    friend ostream& operator << (ostream &s, const Point &p) { return s << '(' << p.x << ", " << p.y << ')'; }
};
inline Point operator + (const Point &p, const Point &q) { return Point(p.x + q.x, p.y + q.y); }
inline Point operator - (const Point &p, const Point &q) { return Point(p.x - q.x, p.y - q.y); }
inline Point operator * (const Point &p, DD a) { return Point(p.x * a, p.y * a); }
inline Point operator * (DD a, const Point &p) { return Point(a * p.x, a * p.y); }
inline Point operator * (const Point &p, const Point &q) { return Point(p.x * q.x - p.y * q.y, p.x * q.y + p.y * q.x); }
inline Point operator / (const Point &p, DD a) { return Point(p.x / a, p.y / a); }
inline Point conj(const Point &p) { return Point(p.x, -p.y); }
inline Point rot(const Point &p, DD ang) { return Point(cos(ang) * p.x - sin(ang) * p.y, sin(ang) * p.x + cos(ang) * p.y); }
inline Point rot90(const Point &p) { return Point(-p.y, p.x); }
inline DD cross(const Point &p, const Point &q) { return p.x * q.y - p.y * q.x; }
inline DD dot(const Point &p, const Point &q) { return p.x * q.x + p.y * q.y; }
inline DD norm(const Point &p) { return dot(p, p); }
inline DD abs(const Point &p) { return sqrt(dot(p, p)); }
inline DD amp(const Point &p) { DD res = atan2(p.y, p.x); if (res < 0) res += PI * 2; return res; }
inline bool eq(const Point &p, const Point &q) { return abs(p - q) < EPS; }
inline bool operator < (const Point &p, const Point &q) { return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y); }
inline bool operator >(const Point &p, const Point &q) { return (abs(p.x - q.x) > EPS ? p.x > q.x : p.y > q.y); }
inline Point operator / (const Point &p, const Point &q) { return p * conj(q) / norm(q); }

/* Line */
struct Line : vector<Point> {
    Line(Point a = Point(0.0, 0.0), Point b = Point(0.0, 0.0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line &l) { return s << '{' << l[0] << ", " << l[1] << '}'; }
};

/* Circle */
struct Circle : Point {
    DD r;
    Circle(Point p = Point(0.0, 0.0), DD r = 0.0) : Point(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle &c) { return s << '(' << c.x << ", " << c.y << ", " << c.r << ')'; }
};



////////////////////////////
// 円や直線の交点
////////////////////////////

Point proj_for_crosspoint(const Point &p, const Line &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
vector<Point> crosspoint(const Line &l, const Line &m) {
    vector<Point> res;
    DD d = cross(m[1] - m[0], l[1] - l[0]);
    if (abs(d) < EPS) return vector<Point>();
    res.push_back(l[0] + (l[1] - l[0]) * cross(m[1] - m[0], m[1] - l[0]) / d);
    return res;
}
vector<Point> crosspoint(const Circle &e, const Circle &f) {
    vector<Point> res;
    DD d = abs(e - f);
    if (d < EPS) return vector<Point>();
    if (d > e.r + f.r + EPS) return vector<Point>();
    if (d < abs(e.r - f.r) - EPS) return vector<Point>();
    DD rcos = (d * d + e.r * e.r - f.r * f.r) / (2.0 * d), rsin;
    if (e.r - abs(rcos) < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point dir = (f - e) / d;
    Point p1 = e + dir * Point(rcos, rsin);
    Point p2 = e + dir * Point(rcos, -rsin);
    res.push_back(p1);
    if (!eq(p1, p2)) res.push_back(p2);
    return res;
}
vector<Point> crosspoint(const Circle &e, const Line &l) {
    vector<Point> res;
    Point p = proj_for_crosspoint(e, l);
    DD rcos = abs(e - p), rsin;
    if (rcos > e.r + EPS) return vector<Point>();
    else if (e.r - rcos < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point dir = (l[1] - l[0]) / abs(l[1] - l[0]);
    Point p1 = p + dir * rsin;
    Point p2 = p - dir * rsin;
    res.push_back(p1);
    if (!eq(p1, p2)) res.push_back(p2);
    return res;
}


// 円と線分の交点
int ccw_for_crosspoint_cs(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
bool isinterPS_crosspoint_cs(const Point &p, const Line &s) {
    return (ccw_for_crosspoint_cs(s[0], s[1], p) == 0);
}
vector<Point> crosspoint_CS(const Circle &e, const Line &s) {
    vector<Point> res;
    Point p = proj_for_crosspoint(e, s);
    DD rcos = abs(e - p), rsin;
    if (rcos > e.r + EPS) return vector<Point>();
    else if (e.r - rcos < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point dir = (s[1] - s[0]) / abs(s[1] - s[0]);
    Point p1 = p - dir * rsin;
    Point p2 = p + dir * rsin;
    if (isinterPS_crosspoint_cs(p1, s)) res.push_back(p1);
    if (isinterPS_crosspoint_cs(p2, s) && !eq(p1, p2)) res.push_back(p2);
    return res;
}


// 円 C と、三角形 ((0, 0), ia, ib) との共通部分の面積
DD calc_common_area(const Circle &c, const Point &ia, const Point &ib) {
    Point a = ia - c, b = ib - c;
    if (abs(a - b) < EPS) return 0;
    auto sub = [&](const Point &x, const Point &y, bool triangle) {
        if (triangle) return cross(x, y) / 2;
        else {
            Point tmp = y * Point(x.x, -x.y);
            DD ang = atan2(tmp.y, tmp.x);
            return c.r * c.r * ang / 2;
        }
    };
    bool isin_a = (abs(a) < c.r + EPS);
    bool isin_b = (abs(b) < c.r + EPS);
    if (isin_a && isin_b) return sub(a, b, true);
    Circle oc(Point(0, 0), c.r);
    Line seg(a, b);
    auto cr = crosspoint_CS(oc, seg);
    if (cr.empty()) return sub(a, b, false);
    auto s = cr[0], t = cr.back();
    return sub(a, s, isin_a) + sub(s, t, true) + sub(t, b, isin_b);
}

DD calc_common_area(const Circle &c, const vector<Point> &pol) {
    DD res = 0;
    int N = pol.size();
    for (int i = 0; i < N; ++i) {
        res += calc_common_area(c, pol[i], pol[(i+1)%N]);
    }
    return res;
}



////////////////////////////
// AOJ 1183 鎖中経路
////////////////////////////

bool isin(const vector<Circle> &cs, const Line &l) {
    vector<Point> vps;
    for (auto c : cs) {
        auto cps = crosspoint_CS(c, l);
        for (auto p : cps) vps.push_back(p);
    }
    sort(vps.begin(), vps.end());
    for (int i = 0; i+1 < vps.size(); ++i) {
        Point p = (vps[i] + vps[i+1]) / 2;
        bool exist = false;
        for (auto c : cs) {
            if (abs(c - p) < c.r + EPS) exist = true;
        }
        if (!exist) return false;
    }
    return true;
}


int main() {
    int N;
    while (cin >> N, N) {
        vector<Circle> cs(N);
        for (int i = 0; i < N; ++i) cin >> cs[i].x >> cs[i].y >> cs[i].r;
        
        vector<Point> ips;
        ips.push_back(Point(cs[0].x, cs[0].y));
        for (int i = 0; i+1 < N; ++i) {
            auto cp = crosspoint(cs[i], cs[i+1]);
            for (auto p : cp) ips.push_back(p);
        }
        ips.push_back(Point(cs[N-1].x, cs[N-1].y));
        
        vector<DD> dp((int)ips.size() + 1, INF);
        dp[0] = 0;
        for (int i = 0; i < ips.size(); ++i) {
            for (int j = i+1; j < ips.size(); ++j) {
                Line l(ips[i], ips[j]);
                if (isin(cs, l)) dp[j] = min(dp[j], dp[i] + abs(l[0] - l[1]));
            }
        }
        cout << fixed << setprecision(10) << dp[(int)ips.size()-1] << endl;
    }
}
