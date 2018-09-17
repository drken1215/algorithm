//
// 3 次元幾何
//
// verified:
//   AOJ 0115 Starship UAZ Advance
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=0115

//   AOJ 1523 Cone Cut
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1523
//


#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


/////////////////////////////////
// 3 次元幾何ライブラリ一式
/////////////////////////////////

using DD = double;
const DD INF = 1LL<<60;      // to be set appropriately
const DD EPS = 1e-10;        // to be set appropriately
const DD PI = acos(-1.0);
DD torad(int deg) {return (DD)(deg) * PI / 180;}
DD todeg(DD ang) {return ang * 180 / PI;}

/* Point */
struct Point3D {
    DD x, y, z;
    Point3D(DD x = 0.0, DD y = 0.0, DD z = 0.0) : x(x), y(y), z(z) {}
    friend ostream& operator << (ostream &s, const Point3D &p) {return s << '(' << p.x << ", " << p.y << ", " << p.z << ')';}
};
Point3D operator + (const Point3D &p, const Point3D &q) {return Point3D(p.x + q.x, p.y + q.y, p.z + q.z);}
Point3D operator - (const Point3D &p, const Point3D &q) {return Point3D(p.x - q.x, p.y - q.y, p.z - q.z);}
Point3D operator * (const Point3D &p, DD a) {return Point3D(p.x * a, p.y * a, p.z * a);}
Point3D operator * (DD a, const Point3D &p) {return Point3D(a * p.x, a * p.y, a * p.z);}
Point3D operator / (const Point3D &p, DD a) {return Point3D(p.x / a, p.y / a), p.z / a;}
Point3D cross (const Point3D &p, const Point3D &q) {
    return Point3D(p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x);
}
DD dot(const Point3D &p, const Point3D &q) {return p.x * q.x + p.y * q.y + p.z * q.z;}
DD norm(const Point3D &p) {return dot(p, p);}
DD abs(const Point3D &p) {return sqrt(dot(p, p));}
bool eq(const Point3D &p, const Point3D &q) {return abs(p - q) < EPS;}
DD area(const Point3D &a, const Point3D &b, const Point3D &c) { return abs(cross(b - a, c - a)) / 2; }

struct Line3D : vector<Point3D> {
    Line3D(const Point3D &a = Point3D(), const Point3D &b = Point3D()) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line3D &l) {return s << '{' << l[0] << ", " << l[1] << '}';}
};

struct Sphere : Point3D {
    DD r;
    Sphere(const Point3D &p = Point3D(), DD r = 0.0) : Point3D(p), r(r) {}
    friend ostream& operator << (ostream &s, const Sphere &c) {return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';}
};

struct Plane : vector<Point3D> {
    Plane(const Point3D &a = Point3D(), const Point3D &b = Point3D(), const Point3D &c = Point3D()) {
        this->push_back(a);
        this->push_back(b);
        this->push_back(c);
    }
    friend ostream& operator << (ostream &s, const Plane &p) {
        return s << '{' << p[0] << ", " << p[1] << ", " << p[2] << '}';
    }
};


Point3D proj(const Point3D &p, const Line3D &l) {
    DD t = dot(p - l[0], l[1] - l[0]) / norm(l[1] - l[0]);
    return l[0] + (l[1] - l[0]) * t;
}

Point3D proj(const Point3D &p, const Plane &pl) {
    Point3D ph = cross(pl[1] - pl[0], pl[2] - pl[0]);
    Point3D pt = proj(p, Line3D(pl[0], pl[0] + ph));
    return p + (pl[0] - pt);
}

Point3D refl(const Point3D &p, const Line3D &l) {
    return p + (proj(p, l) - p) * 2;
}

Point3D refl(const Point3D &p, const Plane &pl) {
    return p + (proj(p, pl) - p) * 2;
}

bool isinterPL(const Point3D &p, const Line3D &l) {
    return (abs(p - proj(p, l)) < EPS);
}
DD distancePL(const Point3D &p, const Line3D &l) {
    return abs(p - proj(p, l));
}
DD distanceLL(const Line3D &l, const Line3D &m) {
    Point3D nv = cross(l[1] - l[0], m[1] - m[0]);
    if (abs(nv) < EPS) return distancePL(l[0], m);
    Point3D p = m[0] - l[0];
    return abs(dot(nv, p)) / abs(nv);
}


vector<Point3D> crosspoint(const Line3D &l, const Plane &pl) {
    vector<Point3D> res;
    Point3D ph = cross(pl[1] - pl[0], pl[2] - pl[0]);
    DD baseLength = dot(l[1] - l[0], ph);
    if (abs(baseLength) < EPS) return vector<Point3D>();
    DD crossLength = dot(pl[0] - l[0], ph);
    DD ratio = crossLength / baseLength;
    res.push_back(l[0] + (l[1] - l[0]) * ratio);
    return res;
}

vector<Point3D> crosspointSPL(const Line3D &s, const Plane &pl) {
    vector<Point3D> res;
    Point3D ph = cross(pl[1] - pl[0], pl[2] - pl[0]);
    DD baseLength = dot(s[1] - s[0], ph);
    if (abs(baseLength) < EPS) return vector<Point3D>();
    DD crossLength = dot(pl[0] - s[0], ph);
    DD ratio = crossLength / baseLength;
    if (ratio < -EPS || ratio > 1.0 + EPS) return vector<Point3D>();
    res.push_back(s[0] + (s[1] - s[0]) * ratio);
    return res;
}



/////////////////////////////////
// AOJ 0115 Starship UAZ Advance
/////////////////////////////////

void solveAOJ0115() {
    Point3D my, en, bar[3];
    cin >> my.x >> my.y >> my.z >> en.x >> en.y >> en.z;
    for (int i = 0; i < 3; ++i) {
        cin >> bar[i].x >> bar[i].y >> bar[i].z;
    }
    Line3D myen(my, en);
    Plane barier(bar[0], bar[1], bar[2]);
    vector<Point3D> cps = crosspointSPL(myen, barier);
    if (cps.empty()) cout << "HIT" << endl;
    else {
        Point3D cp = cps[0];
        if (area(cp, bar[1], bar[2]) + area(cp, bar[2], bar[0]) + area(cp, bar[0], bar[1]) - area(bar[0], bar[1], bar[2]) > EPS)
            cout << "HIT" << endl;
        else
            cout << "MISS" << endl;
    }
}



/////////////////////////////////
// AOJ 1523 Cone Cuts
/////////////////////////////////

// 2 次元
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

struct Line : vector<Point> {
    Line(Point a = Point(0.0, 0.0), Point b = Point(0.0, 0.0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line &l) {return s << '{' << l[0] << ", " << l[1] << '}';}
};

Point proj(const Point &p, const Line &l) {
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

void solveAOJ1523() {
    Point3D X, Y, P;
    DD r;
    cin >> X.x >> X.y >> X.z >> Y.x >> Y.y >> Y.z >> r >> P.x >> P.y >> P.z;
    
    Line3D l(X, Y);
    Point3D PH = proj(P, l);
    Point x(0, abs(X - Y));
    Point y(0, 0);
    Point p(abs(P-PH), abs(PH-Y));
    Point a(r, 0);
    Point b(-r, 0);
        
    vector<Point> vc = crosspoint(Line(p, a), Line(x, b));
    Point c = vc[0];
        
    vector<Point> vd = crosspoint(Line(p, b), Line(x, a));
    Point d = vd[0];
        
    Point m = (c + d)/2;
    Point h = proj(x, Line(c, d));
        
    DD tsr = r * abs(x.y - m.y) / abs(x - y);
    DD sr = sqrt(tsr * tsr - m.x * m.x);
        
    DD tot = PI * r * r * abs(x-y) / 3;
    DD sol = PI * abs(c - d) * sr * abs(x - h) / 6;
    cout << fixed << setprecision(9) << sol << " " << tot-sol << endl;
}



/////////////////////////////////
// main
/////////////////////////////////

int main() {
    solveAOJ0115();
    //solveAOJ1523();
}
