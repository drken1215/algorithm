//
// floor sum
//
// verified
//   Yosupo Library Checker - Sum of Floor of Linear
//     https://judge.yosupo.jp/problem/sum_of_floor_of_linear
//
//   AtCoder ABC 283 Ex - Popcount Sum
//     https://atcoder.jp/contests/abc283/tasks/abc283_h
//
//   AOJ 3217 Cafe au lait (OUPC 2020 I)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=3217
//
//   yukicoder No.2066 Simple Math !
//     https://yukicoder.me/problems/no/2066
//


#include <bits/stdc++.h>
using namespace std;

/*
 sum_{i=0}^{n-1} floor((a * i + b) / m)
 O(log(n + m + a + b))
 
 __int128 can be used for T
 */

// sum_{i=0}^{n-1} floor((a * i + b) / m)
template<class T> T floor_sum(T n, T a, T b, T m) {
    assert(n >= 0 && m >= 1);
    T res = 0;
    if (a < 0) {
        T a2 = (a % m + m) % m;
        res -= n * (n - 1) / 2 * ((a2 - a) / m);
        a = a2;
    }
    if (b < 0) {
        T b2 = (b % m + m) % m;
        res -= n * ((b2 - b) / m);
        b = b2;
    }
    
    while (true) {
        if (a >= m) {
            res += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m) {
            res += n * (b / m);
            b %= m;
        }
        T y_max = a * n + b;
        if (y_max < m) break;
        n = y_max / m;
        b = y_max % m;
        swap(m, a);
    }
    return res;
}

// #lp under (and on) the segment (x1, y1)-(x2, y2)
// not including y = 0, x = x2
template<class T> T num_lattice_points(T x1, T y1, T x2, T y2) {
    T dx = x2 - x1;
    return floor_sum(dx, y2 - y1, dx * y1, dx);
}


// often use
struct i128 {
    // inner value
    __int128 val;
    
    // constructor
    constexpr i128() : val(0) {}
    constexpr i128(long long v) : val(v) {}
    i128(const string &s) : val(0) {
        parse(s);
    }
    void parse(const string &s) {
        val = 0;
        for (auto c : s) {
            if (isdigit(c)) val = val * 10 + (c - '0');
        }
        if (s[0] == '-') val *= -1;
    }
    constexpr __int128 get() const {
        return val;
    }
    constexpr i128 abs() {
        if (val < 0) return -val;
        else return val;
    }
    
    // comparison operators
    constexpr bool operator == (const i128 &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const i128 &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const i128 &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const i128 &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const i128 &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const i128 &r) const {
        return this->val >= r.val;
    }
    
    // arithmetic operators
    constexpr i128& operator += (const i128 &r) {
        val += r.val;
        return *this;
    }
    constexpr i128& operator -= (const i128 &r) {
        val -= r.val;
        return *this;
    }
    constexpr i128& operator *= (const i128 &r) {
        val *= r.val;
        return *this;
    }
    constexpr i128& operator /= (const i128 &r) {
        val /= r.val;
        return *this;
    }
    constexpr i128& operator %= (const i128 &r) {
        val %= r.val;
        return *this;
    }
    constexpr i128 operator + () const {
        return i128(*this);
    }
    constexpr i128 operator - () const {
        return i128(0) - i128(*this);
    }
    constexpr i128 operator + (const i128 &r) const {
        return i128(*this) += r;
    }
    constexpr i128 operator - (const i128 &r) const {
        return i128(*this) -= r;
    }
    constexpr i128 operator * (const i128 &r) const {
        return i128(*this) *= r;
    }
    constexpr i128 operator / (const i128 &r) const {
        return i128(*this) /= r;
    }
    constexpr i128 operator % (const i128 &r) const {
        return i128(*this) %= r;
    }
    
    // bit operators
    constexpr i128 operator >>= (long long r) {
        val >>= r;
        return *this;
    }
    constexpr i128 operator <<= (long long r) {
        val <<= r;
        return *this;
    }
    constexpr i128 operator &= (long long r) {
        val &= r;
        return *this;
    }
    constexpr i128 operator |= (long long r) {
        val |= r;
        return *this;
    }
    constexpr i128 operator << (long long r) const {
        return i128(*this) <<= r;
    }
    constexpr i128 operator >> (long long r) const {
        return i128(*this) >>= r;
    }
    constexpr i128 operator & (long long r) const {
        return i128(*this) &= r;
    }
    constexpr i128 operator | (long long r) const {
        return i128(*this) |= r;
    }
    
    // other operators
    constexpr i128& operator ++ () {
        ++val;
        return *this;
    }
    constexpr i128& operator -- () {
        --val;
        return *this;
    }
    constexpr i128 operator ++ (int) {
        i128 res = *this;
        ++*this;
        return res;
    }
    constexpr i128 operator -- (int) {
        i128 res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, i128 &x) {
        string s;
        is >> s;
        x.parse(s);
        return is;
    }
    friend ostream& operator << (ostream &os, const i128 &x) {
        auto tmp = x.val < 0 ? -x.val : x.val;
        char buffer[128];
        char *d = end(buffer);
        do {
            --d;
            *d = "0123456789"[tmp % 10];
            tmp /= 10;
        } while (tmp != 0);
        if (x.val < 0) {
            --d;
            *d = '-';
        }
        int len = end(buffer) - d;
        if (os.rdbuf()->sputn(d, len) != len) {
            os.setstate(ios_base::badbit);
        }
        return os;
    }
};



//------------------------------//
// Examples
//------------------------------//

// Library Checker
void YosupoSumOfFloorOfLinear() {
    int T;
    cin >> T;
    while (T--) {
        long long N, M, A, B;
        cin >> N >> M >> A >> B;
        cout << floor_sum(N, A, B, M) << endl;
    }
}

// ABC 283 Ex
void ABC_283_Ex() {
    int CASE;
    cin >> CASE;
    while (CASE--) {
        long long N, M, R;
        cin >> N >> M >> R;
        long long range = (N - R) / M + 1;
        
        long long res = 0;
        for (int k = 0; k <= 30; ++k) {
            long long upper = floor_sum(range, M, R + (1LL<<k), 1LL<<(k+1));
            long long lower = floor_sum(range, M, R, 1LL<<(k+1));
            res += upper - lower;
        }
        cout << res << endl;
    }
}

// AOJ 3217
void AOJ_3217() {
    const int MOD = 1000000007;
    struct mint {
        long long val;
        mint(long long v = 0) : val(v % MOD) {
            if (val < 0) val += MOD;
        }
        int getmod() const { return MOD; }
        mint operator - () const {
            return val ? MOD - val : 0;
        }
        mint operator + (const mint& r) const { return mint(*this) += r; }
        mint operator - (const mint& r) const { return mint(*this) -= r; }
        mint operator * (const mint& r) const { return mint(*this) *= r; }
        mint operator / (const mint& r) const { return mint(*this) /= r; }
        mint& operator += (const mint& r) {
            val += r.val;
            if (val >= MOD) val -= MOD;
            return *this;
        }
        mint& operator -= (const mint& r) {
            val -= r.val;
            if (val < 0) val += MOD;
            return *this;
        }
        mint& operator *= (const mint& r) {
            val = val * r.val % MOD;
            return *this;
        }
        mint& operator /= (const mint& r) {
            long long a = r.val, b = MOD, u = 1, v = 0;
            while (b) {
                long long t = a / b;
                a -= t * b, swap(a, b);
                u -= t * v, swap(u, v);
            }
            val = val * u % MOD;
            if (val < 0) val += MOD;
            return *this;
        }
        bool operator == (const mint& r) const {
            return this->val == r.val;
        }
        bool operator != (const mint& r) const {
            return this->val != r.val;
        }
    };

    int N;
    cin >> N;
    vector<long long> X(N), Y(N);
    for (int i = 0; i < N; ++i) cin >> X[i] >> Y[i];
    vector<int> ids(N);
    iota(ids.begin(), ids.end(), 0);
    sort(ids.begin(), ids.end(), [&](int i, int j) {
        return Y[i]*X[j] < Y[j]*X[i];
    });
    mint res = 0;
    long long sy = 0;
    for (auto i : ids) {
        res -= mint(X[i]) * mint(sy);
        res -= num_lattice_points(0LL, 0LL, X[i], Y[i]);
        res += gcd(X[i], Y[i]);
        sy += Y[i];
    }
    reverse(ids.begin(), ids.end());
    sy = 0;
    for (auto i : ids) {
        res += mint(X[i]) * mint(sy);
        res += num_lattice_points(0LL, 0LL, X[i], Y[i]);
        sy += Y[i];
    }
    cout << res.val << endl;
}

// yukicoder No.2066
void Yukicoder_2066() {
    // calc #n (0 <= n <= M) that can be expressed n =  Px + Qy (P, Q is coprime)
    auto calc_num = [&](i128 P, i128 Q, i128 M) {
        i128 mp = M / P;
        i128 N = min(mp + 1, Q);
        i128 a = P, b = M + Q - a * (N - 1);
        return floor_sum(N, a, b, Q) - 1;
    };
    int CASE;
    cin >> CASE;
    while (CASE--) {
        long long P, Q, K;
        cin >> P >> Q >> K;
        long long G = gcd(P, Q);
        P /= G, Q /= G;
        long long low = -1, high = 1LL<<50;
        while (high - low > 1) {
            long long M = (low + high) / 2;
            if (calc_num(P, Q, M) >= K) high = M;
            else low = M;
        }
        cout << high * G << endl;
    }
}


int main() {
    //YosupoSumOfFloorOfLinear();
    //ABC_283_Ex();
    //AOJ_3217();
    Yukicoder_2066();
}
