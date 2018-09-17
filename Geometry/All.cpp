#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


/* Point */
using DD = double;
const DD INF = 1LL<<60;      // to be set appropriately
const DD EPS = 1e-10;        // to be set appropriately
const DD PI = acos(-1.0);
DD torad(int deg) {return (DD)(deg) * PI / 180;}
DD todeg(DD ang) {return ang * 180 / PI;}

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



////////////////////////////
// 点や線分の位置関係
////////////////////////////

// 粗
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)、0：一直線上
int simple_ccw(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    return 0;
}

// 精
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)
// 2：c-a-bの順に一直線上、-2：a-b-cの順に一直線上、0：a-c-bの順に一直線上
int ccw(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < 0) return 2;
    if (norm(b-a) < norm(c-a)) return -2;
    return 0;
}

// 点と三角形の包含関係(辺上については判定していない)
bool is_contain(const Point &p, const Point &a, const Point &b, const Point &c) {
    int r1 = simple_ccw(p, b, c), r2 = simple_ccw(p, c, a), r3 = simple_ccw(p, a, b);
    if (r1 == 1 && r2 == 1 && r3 == 1) return true;
    if (r1 == -1 && r2 == -1 && r3 == -1) return true;
    return false;
}



////////////////////////////
// 特殊な直線, 円を求める
////////////////////////////

// AOJ 2385
// 2点の垂直二等分線（点pを左側に見る向き）
Line bisector(const Point &p, const Point &q) {
    Point c = (p + q) / 2.0L;
    Point v = (q - p) * Point(0.0L, 1.0L);
    v = v / abs(v);
    return Line(c - v, c + v);
}

// AOJ 1039
// 2点の比率a:bのアポロニウスの円
Circle Apporonius(const Point &p, const Point &q, DD a, DD b) {
    if ( abs(a-b) < EPS ) return Circle(Point(0,0),0);
    Point c1 = (p * b + q * a) / (a + b);
    Point c2 = (p * b - q * a) / (b - a);
    Point c = (c1 + c2) / 2;
    DD r = abs(c - c1);
    return Circle(c, r);
}



////////////////////////////
// 円や直線の交差判定, 距離
////////////////////////////

/*
    ccw を用いている
 
    P: point
    L: Line
    S: Segment
 
    distancePL は、「点」と「直線」の距離
    distancePS は、「点」と「線分」の距離
*/

int ccw_for_dis(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < 0) return 2;
    if (norm(b-a) < norm(c-a)) return -2;
    return 0;
}
Point proj(const Point &p, const Line &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
Point refl(const Point &p, const Line &l) {
    return p + (proj(p, l) - p) * 2;
}
bool isinterPL(const Point &p, const Line &l) {
    return (abs(p - proj(p, l)) < EPS);
}
bool isinterPS(const Point &p, const Line &s) {
    return (ccw_for_dis(s[0], s[1], p) == 0);
}
bool isinterLL(const Line &l, const Line &m) {
    return (abs(cross(l[1] - l[0], m[1] - m[0])) > EPS ||
            abs(cross(l[1] - l[0], m[0] - l[0])) < EPS);
}
bool isinterSS(const Line &s, const Line &t) {
    if (eq(s[0], s[1])) return isinterPS(s[0], t);
    if (eq(t[0], t[1])) return isinterPS(t[0], s);
    return (ccw_for_dis(s[0], s[1], t[0]) * ccw_for_dis(s[0], s[1], t[1]) <= 0 &&
            ccw_for_dis(t[0], t[1], s[0]) * ccw_for_dis(t[0], t[1], s[1]) <= 0);
}
DD distancePL(const Point &p, const Line &l) {
    return abs(p - proj(p, l));
}
DD distancePS(const Point &p, const Line &s) {
    Point h = proj(p, s);
    if (isinterPS(h, s)) return abs(p - h);
    return min(abs(p - s[0]), abs(p - s[1]));
}
DD distanceLL(const Line &l, const Line &m) {
    if (isinterLL(l, m)) return 0;
    else return distancePL(m[0], l);
}
DD distanceSS(const Line &s, const Line &t) {
    if (isinterSS(s, t)) return 0;
    else return min(min(distancePS(s[0], t), distancePS(s[1], t)), min(distancePS(t[0], s), distancePS(t[1], s)));
}



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



///////////////////////
// 多角形アルゴリズム
///////////////////////

// 多角形の符号付面積
DD CalcArea(const vector<Point> &pol) {
    DD res = 0.0;
    for (int i = 0; i < pol.size(); ++i) {
        res += cross(pol[i], pol[(i+1)%pol.size()]);
    }
    return res/2.0L;
}

// 凸包 (一直線上の3点を含めない)
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

// 凸包 (一直線上の3点を含める)
vector<Point> ConvexHullCollinearOK(vector<Point> &ps) {
    int n = (int)ps.size();
    vector<Point> res(2*n);
    sort(ps.begin(), ps.end());
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (k >= 2) {
            while (cross(res[k-1] - res[k-2], ps[i] - res[k-2]) < -EPS) {
                --k;
                if (k < 2) break;
            }
        }
        res[k] = ps[i]; ++k;
    }
    int t = k+1;
    for (int i = n-2; i >= 0; --i) {
        if (k >= t) {
            while (cross(res[k-1] - res[k-2], ps[i] - res[k-2]) < -EPS) {
                --k;
                if (k < t) break;
            }
        }
        res[k] = ps[i]; ++k;
    }
    res.resize(k-1);
    return res;
}



///////////////////////
// main
///////////////////////

int main() {
    
}
