//
// 2 円の共通接線
//
// verified:
//   AOJ 2201 Immortal Jewels
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2201
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



///////////////////////
// 接線
///////////////////////

// 点と円
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

// 円と円の共通接線
vector<Line> comtanline(Circle a, Circle b) {
    vector<Line> res;
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
    if (abs(a - b) > a.r + b.r + EPS) {
        Point p = a * b.r + b * a.r;
        p = p * (1.0 / (a.r + b.r));
        vector<Point> bs = tanline(p, a);
        vector<Point> as = tanline(p, b);
        for (int i = 0; i < min(as.size(), bs.size()); ++i) {
            res.push_back(Line(bs[i], as[i]));
        }
    }
    return res;
}




///////////////////////
// ソルバー
///////////////////////

// 距離
Point proj(const Point &p, const Line &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}
DD distancePL(const Point &p, const Line &l) {
    return abs(p - proj(p, l));
}

// カウント
int count(Line l, vector<Circle> vec, vector<DD> d) {
    int res = 0;
    //cout << endl;
    for (int i = 0; i < vec.size(); ++i) {
        DD dis = distancePL(vec[i], l);
        if (dis >= vec[i].r - EPS && dis <= vec[i].r + d[i] + EPS) ++res;
    }
    return res;
}

int main() {
    int N;
    while (cin >> N) {
        if (N == 0) break;
        
        vector<Circle> vec(N);
        vector<DD> d(N);
        for (int i = 0; i < N; ++i) {
            cin >> vec[i].x >> vec[i].y >> vec[i].r >> d[i];
        }
        if (N == 1) cout << 1 << endl;
        else {
            int res = 0;
            for (int i = 0; i < N; ++i) {
                for (int j = i+1; j < N; ++j) {
                    vector<Circle> I(2, vec[i]); I[1].r += d[i];
                    vector<Circle> J(2, vec[j]); J[1].r += d[j];
                    for (int p = 0; p < 2; ++p) {
                        for (int q = 0; q < 2; ++q) {
                            vector<Line> L = comtanline(I[p], J[q]);
                            for (int k = 0; k < L.size(); ++k) {
                                int tmp = count(L[k], vec, d);
                                res = max(res, tmp);
                            }
                        }
                    }
                }
            }
            cout << res << endl;
        }
    }
}
