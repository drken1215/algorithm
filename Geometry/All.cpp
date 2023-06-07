//
// 幾何ライブラリ (二次元)
//

#include <bits/stdc++.h>
using namespace std;


/*/////////////////////////////*/
// 基本要素 (点, 線分, 円)
/*/////////////////////////////*/

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


/*/////////////////////////////*/
// 点や線分の位置関係
/*/////////////////////////////*/

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
// 円や直線の交差判定, 距離
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
vector<Point> convex_hull_colinear(vector<Point> &ps) {
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

// Voronoi-diagram
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
vector<Line> comtanline(Circle a, Circle b) {
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
    sort(ps.begin(), ps.end());
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
// solvers
/*/////////////////////////////*/

int main() {
    
}


