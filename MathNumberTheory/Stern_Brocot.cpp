//
// Stern-Brocot 木上の二分探索
//
// verified
//   Library Checker - Rational Approximation
//     https://judge.yosupo.jp/problem/rational_approximation
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Stern-Brocot Tree
template<class T> struct SternBrocotTree {
    // binary search on Stern-Brocot Tree
    // return {l (= a/b), r (= c/d(} s.t. l: OK, r: NG
    // and a, b, c, d are maximized where a, b, c, d <= lim
    template<class Func> static tuple<T, T, T, T> binary_search(Func check, T lim) {
        assert(check(0, 1));
        assert(!check(1, 0));
        auto rec = [&](auto &&rec, bool which, T &a, T &b, T c, T d) -> void {
            if (a + c > lim || b + d > lim) return;
            if (check(a + c, b + d) == which) {
                a += c, b += d;
                rec(rec, which, a, b, c + c, d + d);
            }
            if (a + c <= lim && b + d <= lim && check(a + c, b + d) == which) a += c, b += d;
        };
        T a = 0, b = 1, c = 1, d = 0;
        while (a + c <= lim && b + d <= lim) {
            rec(rec, true, a, b, c, d);
            rec(rec, false, c, d, a, b);
        }
        return {a, b, c, d};
    }
};


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Rational Approximation
void LibraryCheckerRatilnalApproximation() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    using sbt = SternBrocotTree<long long>;
    int T;
    cin >> T;
    while (T--) {
        long long N, x, y;
        cin >> N >> x >> y;
        auto check = [&](long long a, long long b) -> bool { return a * y <= b * x; };
        auto [a, b, c, d] = sbt::binary_search(check, N);
        if (a * y == b * x) cout << a << ' ' << b << ' ' << a << ' ' << b << '\n';
        else cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
    }
}


int main() {
    LibraryCheckerRatilnalApproximation();
}