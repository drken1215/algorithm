//
// 偏角ソート
//
// verify:
//   Yosupo Library Checker - Sort Points by Argument
//     https://judge.yosupo.jp/problem/sort_points_by_argument
//


#include <bits/stdc++.h>
using namespace std;


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
    constexpr int sign() const {
        if (x >= -EPS && x <= EPS && y >= -EPS && y <= EPS) return 0;
        else if (y < -EPS || (y >= -EPS && y <= EPS && x > EPS)) return -1;
        else return 1;
    }
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
    friend constexpr int sign(const Point &p) {return p.sign();}
    friend constexpr Point rot90(const Point &p) {return p.rot90();}
    friend constexpr Point rot(const Point &p, long long ang) {return p.rot(ang);}
};

// necessary for some functions
template<class DD> constexpr bool operator < (const Point<DD> &p, const Point<DD> &q) {
    return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
}

// Line
template<class DD> struct Line : vector<Point<DD>> {
    Line(Point<DD> a = Point<DD>(0, 0), Point<DD> b = Point<DD>(0, 0)) {
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

// arg sort
// by defining comparison
template<class DD> void arg_sort(vector<Point<DD>> &v) {
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



//------------------------------//
// Examples
//------------------------------//

void Yosupo_Sort_Points_by_Argument() {
    int N;
    cin >> N;
    vector<Point<long long>> vp(N);
    for (int i = 0; i < N; ++i) {
        cin >> vp[i].x >> vp[i].y;
    }
    
    arg_sort(vp);
    for (const auto &p : vp) cout << p.x << " " << p.y << endl;
}

int main() {
    Yosupo_Sort_Points_by_Argument();
}

