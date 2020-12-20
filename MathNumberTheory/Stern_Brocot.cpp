//
// Stern-Brocot 木
//
// verified
//   AOJ 1208 Rational Irrationals
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1208
//


#include <bits/stdc++.h>
using namespace std;

long long P;
long long x, y, u, v;
bool exist_bigger = false;

// consider s/b
void SternBrocot(long long N, long long sl = 0, long long bl = 1, long long sr = 1, long long br = 0) {
    long long s = sl + sr, b = bl + br;
    if (s > N || b > N) return;

    // left (only when bigger case)
    if (s*s > b*b*P) SternBrocot(N, sl, bl, s, b);

    // consider s/b (monotone increasing)
    {
        if (s*s < b*b*P) u = s, v = b;
        else {
            if (!exist_bigger) x = s, y = b, exist_bigger = true;
            return;
        }
    }

    // right
    SternBrocot(N, s, b, sr, br);
}

int main() {
    long long N;
    while (cin >> P >> N, P) {
        exist_bigger = false;
        SternBrocot(N);
        cout << x << "/" << y << " " << u << "/" << v << endl;
    }
}
