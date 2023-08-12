//
// 幾何ライブラリ (二次元)
//
// verify:
//   AOJ Course CGL_1_C - 反時計回り
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_1_C&lang=jp
//
//   AOJ Course CGL_2_D - 距離
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_D&lang=jp
//
//   AOJ Course CGL_4_A - 凸包
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_A&lang=jp
//
//   AOJ Course CGL_7_G - 円の共通接線
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_7_G&lang=ja
//
//   AOJ Course CGL_7_H - 円と多角形の共通部分
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_7_H&lang=ja
//
//   AtCoder ABC 207 D - Congruence Points (Point の複素数としての割り算)
//     https://atcoder.jp/contests/abc207/tasks/abc207_d
//

#include <bits/stdc++.h>
using namespace std;


/*/////////////////////////////*/
// 幾何の基本要素 (点, 線分, 円)
/*/////////////////////////////*/

// basic settings
using DD = long double;
constexpr long double PI = 3.141592653589793238462643383279502884L;
constexpr long double INF = 1LL<<60;  // to be set appropriately
constexpr long double EPS = 1e-10;    // to be set appropriately
long double torad(int deg) {return (long double)(deg) * PI / 180;}
long double todeg(long double ang) {return ang * 180 / PI;}

// Point or Vector
struct Point {
    DD x, y;
    
    // constructor
    constexpr Point() : x(0), y(0) {}
    constexpr Point(DD x, DD y) : x(x), y(y) {}
    
    // various functions
    constexpr Point conj() const {return Point(x, -y);}
    constexpr DD dot(const Point &r) const {return x * r.x + y * r.y;}
    constexpr DD cross(const Point &r) const {return x * r.y - y * r.x;}
    constexpr DD norm() const {return dot(*this);}
    constexpr long double abs() const {return sqrt(norm());}
    constexpr long double amp() const {
        long double res = atan2(y, x);
        if (res < 0) res += PI*2;
        return res;
    }
    constexpr bool eq(const Point &r) const {return (*this - r).abs() <= EPS;}
    constexpr Point rot90() const {return Point(-y, x);}
    constexpr Point rot(long double ang) const {
        return Point(cos(ang) * x - sin(ang) * y, sin(ang) * x + cos(ang) * y);
    }
    
    // arithmetic operators
    constexpr Point operator - () const {return Point(-x, -y);}
    constexpr Point operator + (const Point &r) const {return Point(*this) += r;}
    constexpr Point operator - (const Point &r) const {return Point(*this) -= r;}
    constexpr Point operator * (const Point &r) const {return Point(*this) *= r;}
    constexpr Point operator / (const Point &r) const {return Point(*this) /= r;}
    constexpr Point operator * (DD r) const {return Point(*this) *= r;}
    constexpr Point operator / (DD r) const {return Point(*this) /= r;}
    constexpr Point& operator += (const Point &r) {
        x += r.x, y += r.y;
        return *this;
    }
    constexpr Point& operator -= (const Point &r) {
        x -= r.x, y -= r.y;
        return *this;
    }
    constexpr Point& operator *= (const Point &r) {
        DD tx = x, ty = y;
        x = tx * r.x - ty * r.y;
        y = tx * r.y + ty * r.x;
        return *this;
    }
    constexpr Point& operator /= (const Point &r) {
        return *this *= r.conj() / r.norm();
    }
    constexpr Point& operator *= (DD r) {
        x *= r, y *= r;
        return *this;
    }
    constexpr Point& operator /= (DD r) {
        x /= r, y /= r;
        return *this;
    }

    // friend functions
    friend ostream& operator << (ostream &s, const Point &p) {
        return s << '(' << p.x << ", " << p.y << ')';
    }
    friend constexpr Point conj(const Point &p) {return p.conj();}
    friend constexpr DD dot(const Point &p, const Point &q) {return p.dot(q);}
    friend constexpr DD cross(const Point &p, const Point &q) {return p.cross(q);}
    friend constexpr DD norm(const Point &p) {return p.norm();}
    friend constexpr long double abs(const Point &p) {return p.abs();}
    friend constexpr long double amp(const Point &p) {return p.amp();}
    friend constexpr bool eq(const Point &p, const Point &q) {return p.eq(q);}
    friend constexpr Point rot90(const Point &p) {return p.rot90();}
    friend constexpr Point rot(const Point &p, long long ang) {return p.rot(ang);}
};

// necessary for some functions
constexpr bool operator < (const Point &p, const Point &q) {
    return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
}

// Line
struct Line : vector<Point> {
    Line(Point a = Point(0.0, 0.0), Point b = Point(0.0, 0.0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line &l) {
        return s << '{' << l[0] << ", " << l[1] << '}';
    }
};

// Circle
struct Circle : Point {
    DD r;
    Circle(Point p = Point(0.0, 0.0), DD r = 0.0) : Point(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle &c) {
        return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';
    }
};


/*/////////////////////////////*/
// 点や線分の位置関係
/*/////////////////////////////*/

// 偏角ソート
void amp_sort(vector<Point> &v) {
    sort(v.begin(), v.end(), [&](Point p, Point q) {
        return (abs(amp(p) - amp(q)) > EPS ? amp(p) < amp(q) : norm(p) < norm(q));
    });
}

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
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}

// 点と三角形の包含関係(辺上については判定していない)
bool is_contain(const Point &p, const Point &a, const Point &b, const Point &c) {
    int r1 = simple_ccw(p, b, c), r2 = simple_ccw(p, c, a), r3 = simple_ccw(p, a, b);
    if (r1 == 1 && r2 == 1 && r3 == 1) return true;
    if (r1 == -1 && r2 == -1 && r3 == -1) return true;
    return false;
}


/*/////////////////////////////*/
// 線分の交差判定や距離計算
/*/////////////////////////////*/

int ccw_for_dis(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
Point proj(const Point &p, const Line &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
Point refl(const Point &p, const Line &l) {
    return p + (proj(p, l) - p) * 2;
}
bool is_inter_PL(const Point &p, const Line &l) {
    return (abs(p - proj(p, l)) < EPS);
}
bool is_inter_PS(const Point &p, const Line &s) {
    return (ccw_for_dis(s[0], s[1], p) == 0);
}
bool is_inter_LL(const Line &l, const Line &m) {
    return (abs(cross(l[1] - l[0], m[1] - m[0])) > EPS ||
            abs(cross(l[1] - l[0], m[0] - l[0])) < EPS);
}
bool is_inter_SS(const Line &s, const Line &t) {
    if (eq(s[0], s[1])) return is_inter_PS(s[0], t);
    if (eq(t[0], t[1])) return is_inter_PS(t[0], s);
    return (ccw_for_dis(s[0], s[1], t[0]) * ccw_for_dis(s[0], s[1], t[1]) <= 0 &&
            ccw_for_dis(t[0], t[1], s[0]) * ccw_for_dis(t[0], t[1], s[1]) <= 0);
}
DD distance_PL(const Point &p, const Line &l) {
    return abs(p - proj(p, l));
}
DD distance_PS(const Point &p, const Line &s) {
    Point h = proj(p, s);
    if (is_inter_PS(h, s)) return abs(p - h);
    return min(abs(p - s[0]), abs(p - s[1]));
}
DD distance_LL(const Line &l, const Line &m) {
    if (is_inter_LL(l, m)) return 0;
    else return distance_PL(m[0], l);
}
DD distance_SS(const Line &s, const Line &t) {
    if (is_inter_SS(s, t)) return 0;
    else return min(min(distance_PS(s[0], t), distance_PS(s[1], t)),
                    min(distance_PS(t[0], s), distance_PS(t[1], s)));
}


/*/////////////////////////////*/
// 円や直線の交点
/*/////////////////////////////*/

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


/*/////////////////////////////*/
// 接線
/*/////////////////////////////*/

// tanline
vector<Point> tanline(const Point &p, const Circle &c) {
    vector<Point> res;
    DD d = norm(p - c);
    DD l = d - c.r * c.r;
    if (l < -EPS) return res;
    if (l <= 0.0) l = 0.0;
    Point cq = (p - c) * (c.r * c.r / d);
    Point qs = rot90((p - c) * (c.r * sqrt(l) / d));
    Point s1 = c + cq + qs, s2 = c + cq - qs;
    res.push_back(s1);
    res.push_back(s2);
    return res;
}

// common tanline, a and b must be different!
// Line[0] is tangent point in a
vector<Line> com_tanline(Circle a, Circle b) {
    vector<Line> res;
    // intersect
    if (abs(a - b) > abs(a.r - b.r) + EPS) {
        if (abs(a.r - b.r) < EPS) {
            Point dir = b - a;
            dir = rot90(dir * (a.r / abs(dir)));
            res.push_back(Line(a + dir, b + dir));
            res.push_back(Line(a - dir, b - dir));
        }
        else {
            Point p = a * -b.r + b * a.r;
            p = p * (1.0 / (a.r - b.r));
            vector<Point> bs = tanline(p, a);
            vector<Point> as = tanline(p, b);
            for (int i = 0; i < min(as.size(), bs.size()); ++i) {
                res.push_back(Line(bs[i], as[i]));
            }
        }
    }
    // inscribed
    else if (abs(abs(a - b) - abs(a.r - b.r)) <= EPS) {
        Point dir = b - a;
        if (a.r > b.r) dir = dir * (a.r / abs(dir));
        else dir = dir * (-a.r / abs(dir));
        Point p = a + dir;
        res.push_back(Line(p, p + rot90(dir)));
    }
    // disjoint
    if (abs(a - b) > a.r + b.r + EPS) {
        Point p = a * b.r + b * a.r;
        p = p * (1.0 / (a.r + b.r));
        vector<Point> bs = tanline(p, a);
        vector<Point> as = tanline(p, b);
        for (int i = 0; i < min(as.size(), bs.size()); ++i) {
            res.push_back(Line(bs[i], as[i]));
        }
    }
    // circumscribed
    else if (abs(abs(a - b) - (a.r + b.r)) <= EPS) {
        Point dir = b - a;
        dir = dir * (a.r / abs(dir));
        Point p = a + dir;
        res.push_back(Line(p, p + rot90(dir)));
    }
    return res;
}


/*/////////////////////////////*/
// 多角形アルゴリズム
/*/////////////////////////////*/

// 多角形の符号付面積
DD calc_area(const vector<Point> &pol) {
    DD res = 0.0;
    for (int i = 0; i < pol.size(); ++i) {
        res += cross(pol[i], pol[(i+1)%pol.size()]);
    }
    return res/2.0L;
}

// 点と多角形の包含関係
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


// 凸性判定
int ccw_for_isconvex(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    return 0;
}
bool is_convex(vector<Point> &ps) {
    int n = (int)ps.size();
    for (int i = 0; i < n; ++i) {
        if (ccw_for_isconvex(ps[i], ps[(i+1)%n], ps[(i+2)%n]) == -1) return false;
    }
    return true;
}

// 凸包 (一直線上の3点を含めない)
vector<Point> convex_hull(vector<Point> &ps) {
    int n = (int)ps.size();
    vector<Point> res(2*n);
    auto cmp = [&](Point p, Point q) -> bool {
        return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
    };
    sort(ps.begin(), ps.end(), cmp);
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
vector<Point> convex_hull_colinear(vector<Point> &ps) {
    int n = (int)ps.size();
    vector<Point> res(2*n);
    auto cmp = [&](Point p, Point q) -> bool {
        return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
    };
    sort(ps.begin(), ps.end(), cmp);
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

// Convex Cut
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
vector<Point> convex_cut(const vector<Point> &pol, const Line &l) {
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

// Voronoi Diagram
// pol: outer polygon, ps: points
// find the polygon nearest to ps[ind]
Line bisector(const Point &p, const Point &q) {
    Point c = (p + q) / 2.0L;
    Point v = (q - p) * Point(0.0L, 1.0L);
    v = v / abs(v);
    return Line(c - v, c + v);
}

vector<Point> voronoi(const vector<Point> &pol, const vector<Point> &ps, int ind) {
    vector<Point> res = pol;
    for (int i = 0; i < ps.size(); ++i) {
        if (i == ind) continue;
        Line l = bisector(ps[ind], ps[i]);
        res = convex_cut(res, l);
    }
    return res;
}


/*/////////////////////////////*/
// 面積アルゴリズム
/*/////////////////////////////*/

// 円と円の共通部分の面積
DD calc_common_area(const Circle &p, const Circle &q) {
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

// 交点
int ccw_for_crosspoint_CS(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
bool isinterPS_crosspoint_CS(const Point &p, const Line &s) {
    return (ccw_for_crosspoint_CS(s[0], s[1], p) == 0);
}
Point proj_for_crosspoint_CS(const Point &p, const Line &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
vector<Point> crosspoint_CS(const Circle &e, const Line &s) {
    vector<Point> res;
    Point p = proj_for_crosspoint_CS(e, s);
    DD rcos = abs(e - p), rsin;
    if (rcos > e.r + EPS) return vector<Point>();
    else if (e.r - rcos < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point dir = (s[1] - s[0]) / abs(s[1] - s[0]);
    Point p1 = p - dir * rsin;
    Point p2 = p + dir * rsin;
    if (isinterPS_crosspoint_CS(p1, s)) res.push_back(p1);
    if (isinterPS_crosspoint_CS(p2, s) && !eq(p1, p2)) res.push_back(p2);
    return res;
}

// 原点, 点 x, 点 y とで囲まれる領域の面積 (三角形 ver と扇型 ver)
DD calc_element(const Point &x, const Point &y, DD r, bool triangle) {
    if (triangle) return cross(x, y) / 2;
    else {
        Point tmp = y * Point(x.x, -x.y);
        DD ang = atan2(tmp.y, tmp.x);
        return r * r * ang / 2;
    }
}

// 円 C と、三角形 ((0, 0), ia, ib) との共通部分の面積
DD calc_common_area(const Circle &c, const Point &ia, const Point &ib) {
    Point a = ia - c, b = ib - c;
    if (abs(a - b) < EPS) return 0;
    bool isin_a = (abs(a) < c.r + EPS);
    bool isin_b = (abs(b) < c.r + EPS);
    if (isin_a && isin_b) return calc_element(a, b, c.r, true);
    Circle oc(Point(0, 0), c.r);
    Line seg(a, b);
    auto cr = crosspoint_CS(oc, seg);
    if (cr.empty()) return calc_element(a, b, c.r, false);
    auto s = cr[0], t = cr.back();
    return calc_element(s, t, c.r, true)
        + calc_element(a, s, c.r, isin_a) + calc_element(t, b, c.r, isin_b);
}

// 円 c と多角形 pol の共通部分の面積
DD calc_common_area(const Circle &c, const vector<Point> &pol) {
    DD res = 0;
    int N = pol.size();
    for (int i = 0; i < N; ++i) {
        res += calc_common_area(c, pol[i], pol[(i+1)%N]);
    }
    return res;
}


/*/////////////////////////////*/
// その他
/*/////////////////////////////*/

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
    auto cmp = [&](Point p, Point q) -> bool {
        return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
    };
    sort(ps.begin(), ps.end(), cmp);
    return DivideAndConqur(ps.begin(), n);
}

// 2点の比率a:bのアポロニウスの円 (AOJ 1039)
Circle Apporonius(const Point &p, const Point &q, DD a, DD b) {
    if ( abs(a-b) < EPS ) return Circle(Point(0,0),0);
    Point c1 = (p * b + q * a) / (a + b);
    Point c2 = (p * b - q * a) / (b - a);
    Point c = (c1 + c2) / 2;
    DD r = abs(c - c1);
    return Circle(c, r);
}


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// AOJ CGL_1_C - 反時計回り
void CGL_1_C() {
    Point a, b, c;
    cin >> a.x >> a.y >> b.x >> b.y;
    int N; cin >> N;
    for (int _ = 0; _ < N; ++_) {
        cin >> c.x >> c.y;
        int type = ccw(a, b, c);
        switch (type) {
            case 1: {cout << "COUNTER_CLOCKWISE" << endl; continue; }
            case -1: {cout << "CLOCKWISE" << endl; continue; }
            case 2: {cout << "ONLINE_BACK" << endl; continue; }
            case -2: {cout << "ONLINE_FRONT" << endl; continue; }
            case 0: {cout << "ON_SEGMENT" << endl; continue; }
        }
    }
}

// AOJ CGL_2_D - 距離
void CGL_2_D() {
    int Q;
    cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        Point x1, y1, x2, y2;
        cin >> x1.x >> x1.y >> y1.x >> y1.y >> x2.x >> x2.y >> y2.x >> y2.y;
        Line s(x1, y1), t(x2, y2);
        cout << fixed << setprecision(10) << distance_SS(s, t) << endl;
    }
}

// AOJ CGL_4_A - 凸包
void CGL_4_A() {
    int n;
    cin >> n;
    vector<Point> ps(n);
    for (int i = 0; i < n; ++i) cin >> ps[i].x >> ps[i].y;
    const auto &pol = convex_hull_colinear(ps);
    auto cmp = [&](Point p, Point q) -> bool {
        return (abs(p.y - q.y) > EPS ? p.y < q.y : p.x < q.x);
    };
    Point minv = pol[0];
    int minp = 0;
    for (int i = 0; i < (int)pol.size(); ++i) {
        if (cmp(pol[i], minv)) minv = pol[i], minp = i;
    }
    cout << pol.size() << endl;
    for (int i = 0; i < (int)pol.size(); ++i) {
        int j = (i + minp) % pol.size();
        cout << fixed << setprecision(0) << pol[j].x << " " << pol[j].y << endl;
    }
}

// AOJ CGL_7_G - 共通接線
void CGL_7_G() {
    Circle p, q;
    cin >> p.x >> p.y >> p.r >> q.x >> q.y >> q.r;
    auto l = com_tanline(p, q);
    vector<Point> res;
    for (int i = 0; i < l.size(); ++i) res.push_back(l[i][0]);
    sort(res.begin(), res.end());
    for (int i = 0; i < res.size(); ++i) {
        cout << fixed << setprecision(10) << res[i].x << " " << res[i].y << endl;
    }
}

// AOJ CGL_7_H - 円と多角形の共通部分
void CGL_7_H() {
    int N;
    DD r;
    cin >> N >> r;
    Circle c(Point(0, 0), r);
    vector<Point> pol(N);
    for (int i = 0; i < N; ++i) cin >> pol[i].x >> pol[i].y;
    cout << fixed << setprecision(10) << calc_common_area(c, pol) << endl;
}

// ABC 207 D - Congruence Points (割り算の verify)
void ABC_207_D() {
    int N;
    cin >> N;
    vector<Point> s(N), t(N);
    for (int i = 0; i < N; ++i) cin >> s[i].x >> s[i].y;
    for (int i = 0; i < N; ++i) cin >> t[i].x >> t[i].y;
    
    // 例外処理
    if (N == 1) {
        cout << "Yes" << endl;
        return;
    }

    // s[0] と t[x] が対応するとする
    bool res = false;
    for (int x = 0; x < N; ++x) {
        // 平行移動させて、s[0] と t[x] が原点に来るようにする
        vector<Point> s2(N), t2(N);
        for (int i = 0; i < N; ++i) s2[i] = s[i] - s[0], t2[i] = t[i] - t[x];
        
        // もう 1 点固定する
        for (int y = 0; y < N; ++y) {
            if (x == y) continue;
            if (abs(abs(t2[y]) - abs(s2[1])) > EPS) continue;
            
            // S の各点を回転させる
            vector<Point> s3(N), t3(N);
            Point kaiten = t2[y] / s2[1];
            for (int i = 0; i < N; ++i) {
                s3[i] = s2[i] * kaiten;
                t3[i] = t2[i];
            }
            
            // s3 と t3 が一致するかを判定する
            bool same = true;
            sort(s3.begin(), s3.end());
            sort(t3.begin(), t3.end());
            for (int i = 0; i < N; ++i) if (!eq(s3[i], t3[i])) same = false;
            if (same) res = true;
        }
    }
    cout << (res ? "Yes" : "No") << endl;
}

int main() {
    //CGL_1_C();
    //CGL_2_D();
    //CGL_4_A();
    CGL_7_G();
    //CGL_7_H();
    //ABC_207_D();
}



