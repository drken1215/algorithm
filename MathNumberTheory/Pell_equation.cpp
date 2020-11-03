//
// ペル方程式
//
// verified
//   AOJ 2116 Subdividing a Land
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2116
//
//   Project Euler 140
//     https://projecteuler.net/problem=140
//

/*
  ペル方程式 x^2 - Dy-2 = ±1 の最小解 (x0, y0), (x0-, y0-) を返す
  解なしのときは(-1, -1)を返す
  一般解は、
  1  : (1, 0) から始めて (x0, y0) を次々と mul していく
  -1 : (x0-, y0-) から始めて (x0, y0) を次々と mul していく
*/


#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


// Pell's Equation
pair<long long, long long> Pell(long long D, int c = 1) {
    if (D == 1) return make_pair(-1, -1);
    long long a = D, b = 1;
    while (a > b) {
        a /= 2;
        b *= 2;
    }
    a = b + 1;
    while (a > b) {
        a = b; 
        b = (b + D/b)/2;
    }
    if (a*a == D) return make_pair(-1, -1);
    
    long long a0 = a;
    bool parity = false;
    b = 1;
    long long x2 = 1, x = a, y2 = 0, y = 1, q;
    while (true) {
        b = (D - a*a) / b;
        q = (a0 + a) / b;
        a = q * b - a;
        parity = !parity;
        if (b == 1) break;
        long long tx = x, tx2 = x2, ty = y, ty2 = y2;
        x2 = tx, x = tx * q + tx2; 
        y2 = ty, y = ty * q + ty2;
    }
    long long x0 = x, y0 = y;
    if (!parity) {
        if (c == 1) return make_pair(x, y);
        else return make_pair(-1, -1);
    }
    else if (c == -1) return make_pair(x, y);
    
    long long tx = x, ty = y;
    x = x0 * tx + D * y0 * ty, y = tx * y0 + x0 * ty;
    return make_pair(x, y);
}

pair<long long, long long> PellMul(pair<long long, long long> p, pair<long long, long long> q, long long D) {
    long long f = p.first * q.first + D * p.second * q.second;
    long long s = p.first * q.second + p.second * q.first;
    return make_pair(f, s);
}



///////////////////////////////////////
// solver
///////////////////////////////////////

void AOJ2116() {
    long long N;
    int id = 0;
    while (cin >> N, N) {
        N *= 2;
        long long sq = (long long)(sqrt(N) + 0.5);
        if (sq * sq == N) cout << "Case " << ++id << ": " << sq << " " << 1 << endl;
        else {
            auto res = Pell(N);
            cout << "Case " << ++id << ": " << res.first << " " << res.second << endl;
        }
    }
}

int main() {
    AOJ2116();
}
