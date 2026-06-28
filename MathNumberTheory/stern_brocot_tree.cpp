//
// Stern-Brocot 木
//
// verified
//   Library Checker - Rational Approximation (for 二分探索)
//     https://judge.yosupo.jp/problem/rational_approximation
//
//   Library Checker - Stern–Brocot Tree (for さまざまな木探索)
//     https://judge.yosupo.jp/problem/stern_brocot_tree
//
//   ICPC アジア地区 京都大会 1999 A - Rational Irrationals (AOJ 1208)
//     https://onlinejudge.u-aizu.ac.jp/problems/1208
//
//   AtCoder ABC 294 F - Sugar Water 2 (for 二分探索)
//     https://atcoder.jp/contests/abc294/tasks/abc294_f
//
//   AtCoder ABC 408 G - A/B < p/q < C/D (for 一般的な go_left(), go_right())
//     https://atcoder.jp/contests/abc408/tasks/abc408_g
//
//   AtCode ABC 333 G - Nearest Fraction (for 二分探索)
//     https://atcoder.jp/contests/abc333/tasks/abc333_g
//
//   JOI 2008 春合宿 Day3-2 - 分数 (Fraction)
//     https://atcoder.jp/contests/joisc2008/tasks/joisc2008_fraction
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Stern-Brocot Tree Node
template<class T> struct SBNode {
    // inner values
    T lx, ly, x, y, rx, ry;
    vector<T> path;  // path from root to the node

    // constructors (X, Y must be coprime)
    SBNode() : lx(0), ly(1), x(1), y(1), rx(1), ry(0) {}
    SBNode(T X, T Y) : SBNode() {
        assert(X >= 1 && Y >= 1);
        while (min(X, Y) > 0) {
            if (X > Y) {
                T d = X / Y;
                X -= d * Y;
                go_right(d - (X == 0 ? 1 : 0));
            } else {
                T d = Y / X;
                Y -= d * X;
                go_left(d - (Y == 0 ? 1 : 0));
            }
        }
    }
    SBNode(const vector<T> &Path) : SBNode() {
        for (const T &d : Path) {
            assert(d != 0);
            if (d < 0) go_left(-d);
            if (d > 0) go_right(d);
        }
        assert(path == Path);
    }
    SBNode(const SBNode&) = default;
    SBNode& operator = (const SBNode&) = default;

    // getters
    friend constexpr ostream& operator << (ostream &os, const SBNode &node) {
        return os << "(" << node.x << ", " << node.y << ")"
                  << ", LB: (" << node.lx << ", " << node.ly << ")"
                  << ", UB: (" << node.rx << ", " << node.ry << ")" << '\n';
    }
    constexpr void dump() const {
        cout << "(" << x << ", " << y << ")"
             << ", LB: (" << lx << ", " << ly << ")"
             << ", UB: (" << rx << ", " << ry << ")" << '\n';
        cout << "path: {";
        for (int i = 0; i < (int)path.size(); i++) {
            if (i) cout << ", ";
            cout << path[i];
        }
        cout << "}" << '\n';
    }
    friend constexpr void dump(const SBNode &node) {
        node.dump();
    }
    constexpr T get_depth() const {
        T res = 0;
        for (auto d : path) res += max(d, -d);
        return res;
    }

    // go left (d steps), go right (d steps), go parent (d steps)
    constexpr void go_left(T d = 1) {
        if (d <= 0) return;
        if (path.empty() || path.back() > 0) path.emplace_back(0);
        path.back() -= d;
        rx += lx * d, ry += ly * d;
        x = lx + rx, y = ly + ry;
    }
    constexpr void go_right(T d = 1) {
        if (d <= 0) return;
        if (path.empty() || path.back() < 0) path.emplace_back(0);
        path.back() += d;
        lx += rx * d, ly += ry * d;
        x = lx + rx, y = ly + ry;
    }
    constexpr bool go_parent(T d = 1) {
        if (d <= 0) return true;
        while (d != 0) {
            if (path.empty()) return false;
            T d2 = min(d, max(path.back(), -path.back()));
            if (path.back() > 0) {
                x -= rx * d2, y -= ry * d2;
                lx = x - rx, ly = y - ry;
                path.back() -= d2;
            } else {
                x -= lx * d2, y -= ly * d2;
                rx = x - lx, ry = y - ly;
                path.back() += d2;
            }
            d -= d2;
            if (path.back() == 0) path.pop_back();
            if (d2 == T(0)) break;
        }
        return true;
    }
    friend constexpr void go_left(SBNode &node, T d = 1) { node.go_left(d); }
    friend constexpr void go_right(SBNode &node, T d = 1) { node.go_right(d); }
    friend constexpr bool go_parent(SBNode &node, T d = 1) { return node.go_parent(d); }

    // go left while check(x, y) is False, go right while check(x, y) is True
    template<class Func> constexpr void go_left
    (const Func &check, T lim = numeric_limits<T>::max()/3) {
        assert(check(0, 1));
        assert(!check(1, 0));
        assert(check(lx, ly));
        auto rec = [&](auto &&rec, T dx, T dy, T d = 1) -> void {
            if (rx + dx > lim || ry + dy > lim) return;
            if (!check(rx + dx, ry + dy)) {
                go_left(d);
                rec(rec, dx + dx, dy + dy, d << 1);
            }
            if (rx + dx <= lim && ry + dy <= lim && !check(rx + dx, ry + dy)) {
                go_left(d);
            }
        };
        rec(rec, lx, ly);
    }
    template<class Func> constexpr void go_right
    (const Func &check, T lim = numeric_limits<T>::max()/3) {
        assert(check(0, 1));
        assert(!check(1, 0));
        assert(!check(rx, ry));
        auto rec = [&](auto &&rec, T dx, T dy, T d = 1) -> void {
            if (lx + dx > lim || ly + dy > lim) return;
            if (check(lx + dx, ly + dy)) {
                go_right(d);
                rec(rec, dx + dx, dy + dy, d << 1);
            }
            if (lx + dx <= lim && ly + dy <= lim && check(lx + dx, ly + dy)) {
                go_right(d);
            }
        };
        rec(rec, rx, ry);
    }
    template<class Func> friend constexpr void go_left
    (SBNode &node, const Func &check, T lim = numeric_limits<T>::max()/3) {
        node.go_left(check, lim);
    }
    template<class Func> friend constexpr void go_right
    (SBNode &node, const Func &check, T lim = numeric_limits<T>::max()/3) {
        node.go_right(check, lim);
    }

    // comparison operators
    friend bool operator == (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x == rhs.x && lhs.y == rhs.y;
    }
    friend bool operator != (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x != rhs.x || lhs.y != rhs.y;
    }
    friend bool operator < (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x * rhs.y < rhs.x * lhs.y;
    }
    friend bool operator > (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x * rhs.y > rhs.x * lhs.y;
    }
    friend bool operator <= (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x * rhs.y <= rhs.x * lhs.y;
    }
    friend bool operator >= (const SBNode &lhs, const SBNode &rhs) {
        return lhs.x * rhs.y >= rhs.x * lhs.y;
    }
};

// Stern-Brocot Tree
template<class T> struct SBTree {
    using Node = SBNode<T>;

    // get LCA on Stern-Brocot Tree
    static Node get_lca(const Node &lhs, const Node &rhs) {
        Node v;
        int siz = min<int>(lhs.path.size(), rhs.path.size());
        for (int i = 0; i < siz; i++) {
            T d1 = lhs.path[i], d2 = rhs.path[i];
            if ((d1 < 0) != (d2 < 0)) break;
            if (d1 < 0) v.go_left(min(-d1, -d2));
            if (d1 > 0) v.go_right(min(d1, d2));
            if (d1 != d2) break;
        }
        return v;
    }

    // binary search on Stern-Brocot Tree
    // return {l (= lx/ly), r (= rx/ry)} s.t. l: True, r: False
    // and lx, ly, rx, ry are maximized where lx, ly, rx, ry <= lim
    template<class Func> static Node binary_search(const Func &check, T lim) {
        assert(check(0, 1));
        assert(!check(1, 0));
        Node v;
        while (v.x <= lim && v.y <= lim) {
            // always, check(lx, ly): True, check(rx, ry): False
            v.go_left(check, lim);
            v.go_right(check, lim);
        }
        return v;
    }
};

// not necessary, but often use: rational number
template<class T = long long> struct frac {
    // gcd
    static T gcd(T a, T b) {
        a = max(a, -a), b = max(b, -b);
        while (b) {
            a %= b;
            swap(a, b);
        }
        return a;
    }

    // inner values
    T first, second;

    // constructor
    frac& normalize() {
        if (first == 0 && second != 0) {
            second = 1;
            return *this;
        }
        if (second == 0 && first != 0) {
            first = 1;
            return *this;
        }
        if (second < 0) first = -first, second = -second;
        T d = gcd(max(first, -first), second);
        if (d == 0) first = 0, second = 1;
        else first /= d, second /= d;
        return *this;
    }
    frac(const frac&) = default;
    frac& operator = (const frac&) = default;
    constexpr frac(T f = 0, T s = 1) : first(f), second(s) { 
        normalize(); 
    }
    constexpr frac& operator = (T a) { 
        *this = frac(a, 1); 
        return *this;
    }
    constexpr long double to_double() const {
        assert(second != 0);
        return (long double)(first) / (long double)(second);
    }
    friend constexpr long double to_double(const frac &r) {
        return r.to_double();
    }

    // comparison operators
    constexpr bool operator == (const frac &r) const {
        return this->first == r.first && this->second == r.second;
    }
    constexpr bool operator != (const frac &r) const {
        return this->first != r.first || this->second != r.second;
    }
    constexpr bool operator < (const frac &r) const {
        return this->first * r.second < this->second * r.first;
    }
    constexpr bool operator > (const frac &r) const {
        return this->first * r.second > this->second * r.first;
    }
    constexpr bool operator <= (const frac &r) const {
        return this->first * r.second <= this->second * r.first;
    }
    constexpr bool operator >= (const frac &r) const {
        return this->first * r.second >= this->second * r.first;
    }
    
    // arithmetic operators
    constexpr frac& operator += (const frac &r) {
        this->first = this->first * r.second + this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator -= (const frac &r) {
        this->first = this->first * r.second - this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator *= (const frac &r) {
        this->first *= r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator /= (const frac &r) {
        this->first *= r.second;
        this->second *= r.first;
        this->normalize();
        return *this;
    }
    constexpr frac operator + () const { return frac(*this); }
    constexpr frac operator - () const { return frac(0) - frac(*this); }
    constexpr frac operator + (const frac &r) const { return frac(*this) += r; }
    constexpr frac operator - (const frac &r) const { return frac(*this) -= r; }
    constexpr frac operator * (const frac &r) const { return frac(*this) *= r; }
    constexpr frac operator / (const frac &r) const { return frac(*this) /= r; }
    friend constexpr ostream& operator << (ostream &os, const frac<T> &x) {
        os << x.first; 
        if (x.second != 1) os << "/" << x.second;
        return os;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Rational Approximation
void LibraryCheckerRatilnalApproximation() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    using Node = SBNode<long long>;
    using SBT = SBTree<long long>;
    int T;
    cin >> T;
    while (T--) {
        long long N, x, y;
        cin >> N >> x >> y;
        auto check = [&](long long a, long long b) -> bool { return a * y <= b * x; };
        auto v = SBT::binary_search(check, N);
        if (v.lx * y == v.ly * x) cout << v.lx << ' ' << v.ly << ' ' << v.lx << ' ' << v.ly << '\n';
        else cout << v.lx << ' ' << v.ly << ' ' << v.rx << ' ' << v.ry << '\n';
    }
}

// Library Checker - Stern–Brocot Tree
void LibraryCheckerSternBrocotTree() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    using Node = SBNode<long long>;
    using SBT = SBTree<long long>;
    int T;
    cin >> T;
    while (T--) {
        string qtype;
        cin >> qtype;
        if (qtype == "ENCODE_PATH") {
            long long a, b;
            cin >> a >> b;
            auto res = Node(a, b).path;
            cout << res.size();
            for (auto p : res) {
                if (p > 0) cout << " R " << p;
                else cout << " L " << -p;
            }
            cout << endl;
        } else if (qtype == "DECODE_PATH") {
            long long K;
            cin >> K;
            vector<long long> path(K);
            for (int i = 0; i < K; i++) {
                char c;
                long long n;
                cin >> c >> n;
                if (c == 'L') path[i] = -n;
                else path[i] = n;
            }
            auto res = Node(path);
            cout << res.x << " " << res.y << endl;
        } else if (qtype == "LCA") {
            long long a, b, c, d;
            cin >> a >> b >> c >> d;
            Node l(a, b), r(c, d);
            auto res = SBT::get_lca(l, r);
            cout << res.x << " " << res.y << endl;
        } else if (qtype == "ANCESTOR") {
            long long k, a, b;
            cin >> k >> a >> b;
            Node v(a, b);
            long long d = v.get_depth() - k;
            if (d < 0 || !v.go_parent(d)) cout << -1 << endl;
            else cout << v.x << " " << v.y << endl;
        } else if (qtype == "RANGE") {
            long long a, b;
            cin >> a >> b;
            Node v(a, b);
            cout << v.lx << " " << v.ly << " " << v.rx << " " << v.ry << endl;
        }
    }
}

// ICPC アジア地区 京都大会 1999 A - Rational Irrationals (AOJ 1208)
void AOJ_1208() {
    using SBT = SBTree<long long>;
    long long P, N;
    while (cin >> P >> N, P) {
        auto check = [&](long long a, long long b) -> bool {
            return a*a < b*b*P;
        };
        auto v = SBT::binary_search(check, N);
        cout << v.rx << "/" << v.ry << " " << v.lx << "/" << v.ly << endl;
    }
}

// AtCoder ABC 294 F - Sugar Water 2
void ABC_294_F() {
    using FR = frac<long long>;
    using SBT = SBTree<long long>;
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N), B(N), C(M), D(M);
    for (int i = 0; i < N; i++) cin >> A[i] >> B[i];
    for (int i = 0; i < M; i++) cin >> C[i] >> D[i];

    // x/y 以上の濃度になるものが K 個以上あるか？
    auto check = [&](long long x, long long y) -> bool {
        if (x == 1 && y == 0) return false;
        FR r(x, y);
        vector<FR> P(N), Q(M);
        for (int i = 0; i < N; i++) P[i] = FR(A[i]*100) - r * (A[i]+B[i]);
        for (int i = 0; i < M; i++) Q[i] = FR(C[i]*100) - r * (C[i]+D[i]);
        sort(Q.begin(), Q.end());
        long long num = 0;
        for (int i = 0; i < N; i++) {
            long long tmp = M - (lower_bound(Q.begin(), Q.end(), -P[i]) - Q.begin());
            num += tmp;
        }
        return num >= K;
    };
    auto v = SBT::binary_search(check, 5100000);
    cout << fixed << setprecision(10) << to_double(FR(v.lx, v.ly)) << endl;
}

// AtCoder ABC 408 G - A/B < p/q < C/D
using i128 = __int128_t;
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
         --d;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        --d;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
void ABC_408_G() {
    using Node = SBNode<i128>;
    using SBT = SBTree<i128>;
    int T;
    cin >> T;
    while (T--) {
        i128 A, B, C, D;
        cin >> A >> B >> C >> D;
        auto check_left = [&](i128 x, i128 y) { return A * y >= B * x; };
        auto check_right = [&](i128 x, i128 y) { return x * D < y * C; };
        Node v;
        while (check_left(v.x, v.y) || !check_right(v.x, v.y)) {
            if (check_left(v.x, v.y)) v.go_right(check_left);
            if (!check_right(v.x, v.y)) v.go_left(check_right);
        }
        cout << v.y << endl;
    }
}

// AtCode ABC 333 G - Nearest Fraction
void ABC_333_G() {
    using FR = frac<i128>;
    using SBT = SBTree<i128>;
    string S;
    i128 N, ue = 0, shita = 1;
    cin >> S >> N;
    S = S.substr(2);
    for (auto c : S) {
        ue = ue * 10 + (c - '0');
        shita *= 10;
    }
    FR r(ue, shita);

    auto check = [&](i128 p, i128 q) -> bool { return FR(p, q) <= r; };
    auto v = SBT::binary_search(check, N);
    FR le(v.lx, v.ly), ri(v.rx, v.ry);
    FR dl = max(le - r, r - le), dr = max(ri - r, r - ri);
    if (dl <= dr) cout << v.lx << " " << v.ly << endl;
    else cout << v.rx << " " << v.ry << endl;
}

// JOI 2008 春合宿 Day3-2 - 分数 (Fraction)
void JOI_2008_Day3_2() {
    using Node = SBNode<long long>;
    long long M, K;
    cin >> M >> K;

    // 分母が M 以下の規約分数を小さい順に走査していく
    Node res, v;
    auto rec = [&](auto &&rec, Node &v) -> void {
        if (K <= 0) return;
        v.go_left();
        if (v.y <= M) rec(rec, v);
        v.go_parent();
        K--;
        if (K == 0) {
            res = v;
            return;
        }
        v.go_right();
        if (v.y <= M) rec(rec, v);
        v.go_parent();
    };
    rec(rec, v);
    if (res.x < res.y) cout << res.x << " " << res.y << endl;
    else cout << -1 << endl;
}


int main() {
    //LibraryCheckerRatilnalApproximation();
    //LibraryCheckerSternBrocotTree();
    //AOJ_1208();
    //ABC_294_F();
    //ABC_408_G();
    ABC_333_G();
    //JOI_2008_Day3_2();
}