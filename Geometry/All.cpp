//
// 幾何ライブラリ (二次元)
//
// verify:
//   Yosupo Library Checker - Sort Points by Argument
//     https://judge.yosupo.jp/problem/sort_points_by_argument
//
//   AOJ Course CGL_1_C - 反時計回り
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_1_C&lang=jp
//
//   AOJ Course CGL_2_C - 交点
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_C&lang=ja
//
//   AOJ Course CGL_2_D - 距離
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_D&lang=jp
//
//   AOJ Course CGL_4_A - 凸包
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_A&lang=jp
//
//   AOJ Course CGL_5_A - 最近点対
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_5_A&lang=ja
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
long double EPS = 1e-10;  // to be set appropriately
constexpr long double PI = 3.141592653589793238462643383279502884L;
long double torad(long double deg) {return (long double)(deg) * PI / 180;}
long double todeg(long double ang) {return ang * 180 / PI;}

// Point or Vector
template<class DD> struct Point {
    // inner value
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
    constexpr long double arg() const {
        if (this->eq(Point(0, 0))) return 0L;
        else if (x < -EPS && this->eq(Point(x, 0))) return PI;
        else return atan2((long double)(y), (long double)(x));
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
    friend constexpr long double arg(const Point &p) {return p.arg();}
    friend constexpr bool eq(const Point &p, const Point &q) {return p.eq(q);}
    friend constexpr Point rot90(const Point &p) {return p.rot90();}
    friend constexpr Point rot(const Point &p, long long ang) {return p.rot(ang);}
};

// necessary for some functions
template<class DD> constexpr bool operator < (const Point<DD> &p, const Point<DD> &q) {
    return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
}

// Line
template<class DD> struct Line : vector<Point<DD>> {
    Line(Point<DD> a = Point(0, 0), Point<DD> b = Point(0, 0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line<DD> &l) {
        return s << '{' << l[0] << ", " << l[1] << '}';
    }
};

// Circle
template<class DD> struct Circle : Point<DD> {
    DD r;
    Circle(Point<DD> p = Point<DD>(0, 0), DD r = 0) : Point<DD>(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle<DD> &c) {
        return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';
    }
};


/*/////////////////////////////*/
// 点や線分の位置関係
/*/////////////////////////////*/

// arg sort
// by defining comparison
template<class DD> void arg_sort(vector<Point<DD>> &v) {
    auto sign = [&](const Point<DD> &p) -> int {
        if (abs(p.x) <= EPS && abs(p.y) <= EPS) return 0;
        else if (p.y < -EPS || (abs(p.y) <= EPS && p.x > EPS)) return -1;
        else return 1;
    };
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
        if (sign(p) != sign(q)) return sign(p) < sign(q);
        return (abs(cross(p, q)) > EPS ? cross(p, q) > EPS : norm(p) < norm(q));
    };
    sort(v.begin(), v.end(), cmp);
}
// by calculating arg directly
template<class DD> void arg_sort_direct(vector<Point<DD>> &v) {
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
        return (abs(arg(p) - arg(q)) > EPS ? arg(p) < arg(q) : norm(p) < norm(q));
    };
    sort(v.begin(), v.end(), cmp);
}

// 粗
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)、0：一直線上
template<class DD> int simple_ccw
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    return 0;
}

// 精
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)
// 2：c-a-bの順に一直線上、-2：a-b-cの順に一直線上、0：a-c-bの順に一直線上
template<class DD> int ccw
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}

// 点と三角形の包含関係(辺上については判定していない)
template<class DD> bool is_contain
 (const Point<DD> &p, const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    int r1 = simple_ccw(p, b, c), r2 = simple_ccw(p, c, a), r3 = simple_ccw(p, a, b);
    if (r1 == 1 && r2 == 1 && r3 == 1) return true;
    if (r1 == -1 && r2 == -1 && r3 == -1) return true;
    return false;
}


/*/////////////////////////////*/
// 線分の交差判定や距離計算
/*/////////////////////////////*/

template<class DD> int ccw_for_dis
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
template<class DD> Point<DD> proj(const Point<DD> &p, const Line<DD> &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
template<class DD> Point<DD> refl(const Point<DD> &p, const Line<DD> &l) {
    return p + (proj(p, l) - p) * 2;
}
template<class DD> bool is_inter_PL(const Point<DD> &p, const Line<DD> &l) {
    return (abs(p - proj(p, l)) < EPS);
}
template<class DD> bool is_inter_PS(const Point<DD> &p, const Line<DD> &s) {
    return (ccw_for_dis(s[0], s[1], p) == 0);
}
template<class DD> bool is_inter_LL(const Line<DD> &l, const Line<DD> &m) {
    return (abs(cross(l[1] - l[0], m[1] - m[0])) > EPS ||
            abs(cross(l[1] - l[0], m[0] - l[0])) < EPS);
}
template<class DD> bool is_inter_SS(const Line<DD> &s, const Line<DD> &t) {
    if (eq(s[0], s[1])) return is_inter_PS(s[0], t);
    if (eq(t[0], t[1])) return is_inter_PS(t[0], s);
    return (ccw_for_dis(s[0], s[1], t[0]) * ccw_for_dis(s[0], s[1], t[1]) <= 0 &&
            ccw_for_dis(t[0], t[1], s[0]) * ccw_for_dis(t[0], t[1], s[1]) <= 0);
}
template<class DD> DD distance_PL(const Point<DD> &p, const Line<DD> &l) {
    return abs(p - proj(p, l));
}
template<class DD> DD distance_PS(const Point<DD> &p, const Line<DD> &s) {
    Point h = proj(p, s);
    if (is_inter_PS(h, s)) return abs(p - h);
    return min(abs(p - s[0]), abs(p - s[1]));
}
template<class DD> DD distance_LL(const Line<DD> &l, const Line<DD> &m) {
    if (is_inter_LL(l, m)) return 0;
    else return distance_PL(m[0], l);
}
template<class DD> DD distance_SS(const Line<DD> &s, const Line<DD> &t) {
    if (is_inter_SS(s, t)) return 0;
    else return min(min(distance_PS(s[0], t), distance_PS(s[1], t)),
                    min(distance_PS(t[0], s), distance_PS(t[1], s)));
}


/*/////////////////////////////*/
// 円や直線の交点
/*/////////////////////////////*/

template<class DD> Point<DD> proj_for_crosspoint(const Point<DD> &p, const Line<DD> &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
template<class DD> vector<Point<DD>> crosspoint(const Line<DD> &l, const Line<DD> &m) {
    vector<Point<DD>> res;
    DD d = cross(m[1] - m[0], l[1] - l[0]);
    if (abs(d) < EPS) return vector<Point<DD>>();
    res.push_back(l[0] + (l[1] - l[0]) * cross(m[1] - m[0], m[1] - l[0]) / d);
    return res;
}
template<class DD> vector<Point<DD>> crosspoint_SS(const Line<DD> &l, const Line<DD> &m) {
    if (is_inter_SS(l, m)) return crosspoint(l, m);
    else return vector<Point<DD>>();
}
template<class DD> vector<Point<DD>> crosspoint(const Circle<DD> &e, const Circle<DD> &f) {
    vector<Point<DD>> res;
    DD d = abs(e - f);
    if (d < EPS) return vector<Point<DD>>();
    if (d > e.r + f.r + EPS) return vector<Point<DD>>();
    if (d < abs(e.r - f.r) - EPS) return vector<Point<DD>>();
    DD rcos = (d * d + e.r * e.r - f.r * f.r) / (2.0 * d), rsin;
    if (e.r - abs(rcos) < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point<DD> dir = (f - e) / d;
    Point<DD> p1 = e + dir * Point(rcos, rsin);
    Point<DD> p2 = e + dir * Point(rcos, -rsin);
    res.push_back(p1);
    if (!eq(p1, p2)) res.push_back(p2);
    return res;
}
template<class DD> vector<Point<DD>> crosspoint(const Circle<DD> &e, const Line<DD> &l) {
    vector<Point<DD>> res;
    Point<DD> p = proj_for_crosspoint(e, l);
    DD rcos = abs(e - p), rsin;
    if (rcos > e.r + EPS) return vector<Point<DD>>();
    else if (e.r - rcos < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point<DD> dir = (l[1] - l[0]) / abs(l[1] - l[0]);
    Point<DD> p1 = p + dir * rsin;
    Point<DD> p2 = p - dir * rsin;
    res.push_back(p1);
    if (!eq(p1, p2)) res.push_back(p2);
    return res;
}


/*/////////////////////////////*/
// 接線
/*/////////////////////////////*/

// tanline
template<class DD> vector<Point<DD>> tanline(const Point<DD> &p, const Circle<DD> &c) {
    vector<Point<DD>> res;
    DD d = norm(p - c);
    DD l = d - c.r * c.r;
    if (l < -EPS) return res;
    if (l <= 0.0) l = 0.0;
    Point<DD> cq = (p - c) * (c.r * c.r / d);
    Point<DD> qs = rot90((p - c) * (c.r * sqrt(l) / d));
    Point<DD> s1 = c + cq + qs, s2 = c + cq - qs;
    res.push_back(s1);
    res.push_back(s2);
    return res;
}

// common tanline, a and b must be different!
// Line[0] is tangent point in a
template<class DD> vector<Line<DD>> com_tanline(const Circle<DD> &a, const Circle<DD> &b) {
    vector<Line<DD>> res;
    // intersect
    if (abs(a - b) > abs(a.r - b.r) + EPS) {
        if (abs(a.r - b.r) < EPS) {
            Point<DD> dir = b - a;
            dir = rot90(dir * (a.r / abs(dir)));
            res.push_back(Line(a + dir, b + dir));
            res.push_back(Line(a - dir, b - dir));
        }
        else {
            Point<DD> p = a * -b.r + b * a.r;
            p = p * (1.0 / (a.r - b.r));
            vector<Point<DD>> bs = tanline(p, a);
            vector<Point<DD>> as = tanline(p, b);
            for (int i = 0; i < min(as.size(), bs.size()); ++i) {
                res.push_back(Line(bs[i], as[i]));
            }
        }
    }
    // inscribed
    else if (abs(abs(a - b) - abs(a.r - b.r)) <= EPS) {
        Point<DD> dir = b - a;
        if (a.r > b.r) dir = dir * (a.r / abs(dir));
        else dir = dir * (-a.r / abs(dir));
        Point<DD> p = a + dir;
        res.push_back(Line(p, p + rot90(dir)));
    }
    // disjoint
    if (abs(a - b) > a.r + b.r + EPS) {
        Point<DD> p = a * b.r + b * a.r;
        p = p * (1.0 / (a.r + b.r));
        vector<Point<DD>> bs = tanline(p, a);
        vector<Point<DD>> as = tanline(p, b);
        for (int i = 0; i < min(as.size(), bs.size()); ++i) {
            res.push_back(Line(bs[i], as[i]));
        }
    }
    // circumscribed
    else if (abs(abs(a - b) - (a.r + b.r)) <= EPS) {
        Point<DD> dir = b - a;
        dir = dir * (a.r / abs(dir));
        Point<DD> p = a + dir;
        res.push_back(Line(p, p + rot90(dir)));
    }
    return res;
}


/*/////////////////////////////*/
// 多角形アルゴリズム
/*/////////////////////////////*/

// 多角形の符号付面積
template<class DD> DD calc_area(const vector<Point<DD>> &pol) {
    DD res = 0.0;
    for (int i = 0; i < pol.size(); ++i) {
        res += cross(pol[i], pol[(i+1)%pol.size()]);
    }
    return res/2.0L;
}

// 点と多角形の包含関係
// 2: in, 1: on, 0: out
template<class DD> int is_contain(const vector<Point<DD>> &pol, const Point<DD> &p) {
    int n = (int)pol.size();
    int isin = 0;
    for (int i = 0; i < n; ++i) {
        Point<DD> a = pol[i] - p, b = pol[(i+1)%n] - p;
        if (a.y > b.y) swap(a, b);
        if (a.y <= 0 && b.y > 0) if (cross(a, b) < 0) isin = 1-isin;
        if (cross(a, b) == 0 && dot(a, b) <= 0) return 1;
    }
    if (isin) return 2;
    else return 0;
}


// 凸性判定
template<class DD> int ccw_for_isconvex
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    return 0;
}
template<class DD> bool is_convex(const vector<Point<DD>> &ps) {
    int n = (int)ps.size();
    for (int i = 0; i < n; ++i) {
        if (ccw_for_isconvex(ps[i], ps[(i+1)%n], ps[(i+2)%n]) == -1) return false;
    }
    return true;
}

// 凸包 (一直線上の3点を含めない)
template<class DD> vector<Point<DD>> convex_hull(vector<Point<DD>> &ps) {
    int n = (int)ps.size();
    vector<Point<DD>> res(2*n);
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
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
template<class DD> vector<Point<DD>> convex_hull_colinear(vector<Point<DD>> &ps) {
    int n = (int)ps.size();
    vector<Point<DD>> res(2*n);
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
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
template<class DD> int ccw_for_convexcut
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
template<class DD> vector<Point<DD>> crosspoint_for_convexcut
 (const Line<DD> &l, const Line<DD> &m) {
    vector<Point<DD>> res;
    DD d = cross(m[1] - m[0], l[1] - l[0]);
    if (abs(d) < EPS) return vector<Point<DD>>();
    res.push_back(l[0] + (l[1] - l[0]) * cross(m[1] - m[0], m[1] - l[0]) / d);
    return res;
}
template<class DD> vector<Point<DD>> convex_cut
 (const vector<Point<DD>> &pol, const Line<DD> &l) {
    vector<Point<DD>> res;
    for (int i = 0; i < pol.size(); ++i) {
        Point<DD> p = pol[i], q = pol[(i+1)%pol.size()];
        if (ccw_for_convexcut(l[0], l[1], p) != -1) {
            if (res.size() == 0) res.push_back(p);
            else if (!eq(p, res[res.size()-1])) res.push_back(p);
        }
        if (ccw_for_convexcut(l[0], l[1], p) * ccw_for_convexcut(l[0], l[1], q) < 0) {
            vector<Point<DD>> temp = crosspoint_for_convexcut(Line(p, q), l);
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
template<class DD> Line<DD> bisector(const Point<DD> &p, const Point<DD> &q) {
    Point<DD> c = (p + q) / 2.0L;
    Point<DD> v = (q - p) * Point(0.0L, 1.0L);
    v = v / abs(v);
    return Line(c - v, c + v);
}

template<class DD> vector<Point<DD>> voronoi
 (const vector<Point<DD>> &pol, const vector<Point<DD>> &ps, int ind) {
    vector<Point<DD>> res = pol;
    for (int i = 0; i < ps.size(); ++i) {
        if (i == ind) continue;
        Line<DD> l = bisector(ps[ind], ps[i]);
        res = convex_cut(res, l);
    }
    return res;
}


/*/////////////////////////////*/
// 面積アルゴリズム
/*/////////////////////////////*/

// 円と円の共通部分の面積
template<class DD> DD calc_common_area(const Circle<DD> &p, const Circle<DD> &q) {
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
template<class DD> int ccw_for_crosspoint_CS
 (const Point<DD> &a, const Point<DD> &b, const Point<DD> &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
template<class DD> bool isinterPS_crosspoint_CS(const Point<DD> &p, const Line<DD> &s) {
    return (ccw_for_crosspoint_CS(s[0], s[1], p) == 0);
}
template<class DD> Point<DD> proj_for_crosspoint_CS(const Point<DD> &p, const Line<DD> &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
template<class DD> vector<Point<DD>> crosspoint_CS(const Circle<DD> &e, const Line<DD> &s) {
    vector<Point<DD>> res;
    Point<DD> p = proj_for_crosspoint_CS(e, s);
    DD rcos = abs(e - p), rsin;
    if (rcos > e.r + EPS) return vector<Point<DD>>();
    else if (e.r - rcos < EPS) rsin = 0;
    else rsin = sqrt(e.r * e.r - rcos * rcos);
    Point<DD> dir = (s[1] - s[0]) / abs(s[1] - s[0]);
    Point<DD> p1 = p - dir * rsin;
    Point<DD> p2 = p + dir * rsin;
    if (isinterPS_crosspoint_CS(p1, s)) res.push_back(p1);
    if (isinterPS_crosspoint_CS(p2, s) && !eq(p1, p2)) res.push_back(p2);
    return res;
}

// 原点, 点 x, 点 y とで囲まれる領域の面積 (三角形 ver と扇型 ver)
template<class DD> DD calc_element
 (const Point<DD> &x, const Point<DD> &y, DD r, bool triangle = true) {
    if (triangle) return cross(x, y) / 2;
    else {
        Point<DD> tmp = y * Point(x.x, -x.y);
        DD ang = atan2(tmp.y, tmp.x);
        return r * r * ang / 2;
    }
}

// 円 C と、三角形 ((0, 0), ia, ib) との共通部分の面積
template<class DD> DD calc_common_area
 (const Circle<DD> &c, const Point<DD> &ia, const Point<DD> &ib) {
    Point<DD> a = ia - c, b = ib - c;
    if (abs(a - b) < EPS) return 0;
    bool isin_a = (abs(a) < c.r + EPS);
    bool isin_b = (abs(b) < c.r + EPS);
    if (isin_a && isin_b) return calc_element(a, b, c.r, true);
    Circle<DD> oc(Point<DD>(0, 0), c.r);
    Line<DD> seg(a, b);
    auto cr = crosspoint_CS(oc, seg);
    if (cr.empty()) return calc_element(a, b, c.r, false);
    auto s = cr[0], t = cr.back();
    return calc_element(s, t, c.r, true)
        + calc_element(a, s, c.r, isin_a) + calc_element(t, b, c.r, isin_b);
}

// 円 c と多角形 pol の共通部分の面積
template<class DD> DD calc_common_area(const Circle<DD> &c, const vector<Point<DD>> &pol) {
    DD res = 0;
    int n = (int)pol.size();
    for (int i = 0; i < n; ++i) {
        res += calc_common_area(c, pol[i], pol[(i+1)%n]);
    }
    return res;
}


/*/////////////////////////////*/
// その他
/*/////////////////////////////*/

// 最近点対
template<class DD> DD Cloest(vector<Point<DD>> &ps) {
    typedef typename vector<Point<DD>>::iterator Iterator;
    auto dac = [&](auto self, const Iterator &it, int n) -> DD {
        if (n <= 1) return numeric_limits<DD>::max();
        int m = n/2;
        DD x = it[m].x;
        DD d = min(self(self, it, m), self(self, it+m, n-m));
        inplace_merge(it, it+m, it+n, [&](const Point<DD> &a, const Point<DD> &b){return a.y < b.y;});
        
        vector<Point<DD>> vec;
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
    };
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
        return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
    };
    sort(ps.begin(), ps.end(), cmp);
    return dac(dac, ps.begin(), (int)ps.size());
}

// 2点の比率a:bのアポロニウスの円 (AOJ 1039)
template<class DD> Circle<DD> Apporonius(const Point<DD> &p, const Point<DD> &q, DD a, DD b) {
    if (abs(a-b) < EPS) return Circle<DD>(Point<DD>(0, 0), 0);
    Point<DD> c1 = (p * b + q * a) / (a + b);
    Point<DD> c2 = (p * b - q * a) / (b - a);
    Point<DD> c = (c1 + c2) / 2;
    DD r = abs(c - c1);
    return Circle<DD>(c, r);
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// arg sort
void Yosupo_Sort_Points_by_Argument_1() {
    auto prev_EPS = EPS;
    EPS = 1e-18;  // 10^9 サイズの座標間の角度は 10^-18 オーダーになる
    int N;
    cin >> N;
    vector<Point<long long>> vp(N);
    for (int i = 0; i < N; ++i) cin >> vp[i].x >> vp[i].y;
    
    arg_sort_direct(vp);
    for (const auto &p : vp) cout << p.x << " " << p.y << endl;
    EPS = prev_EPS;
}
void Yosupo_Sort_Points_by_Argument_2() {
    int N;
    cin >> N;
    vector<Point<long long>> vp(N);
    for (int i = 0; i < N; ++i) cin >> vp[i].x >> vp[i].y;
    
    arg_sort(vp);
    for (const auto &p : vp) cout << p.x << " " << p.y << endl;
}

// AOJ CGL_1_C - 反時計回り
void CGL_1_C() {
    using DD = long double;
    Point<DD> a, b, c;
    cin >> a.x >> a.y >> b.x >> b.y;
    int N;
    cin >> N;
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

// AOJ CGL_2_C - 交点
void CGL_2_C() {
    using DD = long double;
    int N;
    cin >> N;
    while (N--) {
        Point<DD> a, b, c, d;
        cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;
        Line<DD> s(a, b), t(c, d);
        const auto &res = crosspoint_SS(s, t);
        cout << fixed << setprecision(10) << res[0].x << " " << res[0].y << endl;
    }
}

// AOJ CGL_2_D - 距離
void CGL_2_D() {
    using DD = long double;
    int Q;
    cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        Point<DD> x1, y1, x2, y2;
        cin >> x1.x >> x1.y >> y1.x >> y1.y >> x2.x >> x2.y >> y2.x >> y2.y;
        Line s(x1, y1), t(x2, y2);
        cout << fixed << setprecision(10) << distance_SS(s, t) << endl;
    }
}

// AOJ CGL_4_A - 凸包
void CGL_4_A() {
    using DD = long double;
    int n;
    cin >> n;
    vector<Point<DD>> ps(n);
    for (int i = 0; i < n; ++i) cin >> ps[i].x >> ps[i].y;
    const auto &pol = convex_hull_colinear(ps);
    auto cmp = [&](const Point<DD> &p, const Point<DD> &q) -> bool {
        return (abs(p.y - q.y) > EPS ? p.y < q.y : p.x < q.x);
    };
    Point<DD> minv = pol[0];
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

// CGL_5_A - 最近点対
void CGL_5_A() {
    using DD = long double;
    int N;
    cin >> N;
    vector<Point<DD>> ps(N);
    for (int i = 0; i < N; ++i) cin >> ps[i].x >> ps[i].y;
    cout << fixed << setprecision(10) << Cloest(ps) << endl;
}

// AOJ CGL_7_G - 共通接線
void CGL_7_G() {
    using DD = long double;
    Circle<DD> p, q;
    cin >> p.x >> p.y >> p.r >> q.x >> q.y >> q.r;
    auto l = com_tanline(p, q);
    vector<Point<DD>> res;
    for (int i = 0; i < l.size(); ++i) res.push_back(l[i][0]);
    sort(res.begin(), res.end());
    for (int i = 0; i < res.size(); ++i) {
        cout << fixed << setprecision(10) << res[i].x << " " << res[i].y << endl;
    }
}

// AOJ CGL_7_H - 円と多角形の共通部分
void CGL_7_H() {
    using DD = long double;
    int N;
    DD r;
    cin >> N >> r;
    Circle<DD> c(Point<DD>(0, 0), r);
    vector<Point<DD>> pol(N);
    for (int i = 0; i < N; ++i) cin >> pol[i].x >> pol[i].y;
    cout << fixed << setprecision(10) << calc_common_area(c, pol) << endl;
}

// ABC 207 D - Congruence Points (割り算の verify)
void ABC_207_D() {
    using DD = long double;
    int N;
    cin >> N;
    vector<Point<DD>> s(N), t(N);
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
        vector<Point<DD>> s2(N), t2(N);
        for (int i = 0; i < N; ++i) s2[i] = s[i] - s[0], t2[i] = t[i] - t[x];
        
        // もう 1 点固定する
        for (int y = 0; y < N; ++y) {
            if (x == y) continue;
            if (abs(abs(t2[y]) - abs(s2[1])) > EPS) continue;
            
            // S の各点を回転させる
            vector<Point<DD>> s3(N), t3(N);
            Point<DD> kaiten = t2[y] / s2[1];
            for (int i = 0; i < N; ++i) {
                s3[i] = s2[i] * kaiten;
                t3[i] = t2[i];
            }
            
            // s3 と t3 が一致するかを判定する
            bool same = true;
            arg_sort(s3), arg_sort(t3);;
            for (int i = 0; i < N; ++i) if (!eq(s3[i], t3[i])) same = false;
            if (same) res = true;
        }
    }
    cout << (res ? "Yes" : "No") << endl;
}


int main() {
    //Yosupo_Sort_Points_by_Argument_1();
    //Yosupo_Sort_Points_by_Argument_2();
    //CGL_1_C();
    //CGL_2_C();
    //CGL_2_D();
    //CGL_4_A();
    //CGL_5_A();
    //CGL_7_G();
    //CGL_7_H();
    ABC_207_D();
}

