//
// 128 ビット整数 __int128 型のラッパー
//
// verified
//   Yosupo Library Checker - Sum of Floor of Linear
//     https://judge.yosupo.jp/problem/sum_of_floor_of_linear
//


#include <bits/stdc++.h>
using namespace std;


// int 128
using i128 = __int128;
constexpr i128 to_integer(const string &s) {
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


//------------------------------//
// Examples
//------------------------------//

void mini_test() {
    i128 a = 99;
    i128 b = to_integer("993421434234324234234432421234324");
    
    cout << a + b << endl;
    cout << a * b << endl;
    cout << a % b << endl;
    cout << b * 2 << endl;
}

// Library Checker
// sum_{i=0}^{n-1} floor((a * i + b) / m)
template<class T> T floor_sum(T n, T a, T b, T m) {
    if (n == 0) return 0;
    T res = 0;
    if (a >= m) {
        res += n * (n - 1) * (a / m) / 2;
        a %= m;
    }
    if (b >= m) {
        res += n * (b / m);
        b %= m;
    }
    if (a == 0) return res;
    T ymax = (a * n + b) / m, xmax = ymax * m - b;
    if (ymax == 0) return res;
    res += (n - (xmax + a - 1) / a) * ymax;
    res += floor_sum(ymax, m, (a - xmax % a) % a, a);
    return res;
}

void YosupoSumOfFloorOfLinear() {
    int T;
    cin >> T;
    while (T--) {
        i128 N, M, A, B;
        cin >> N >> M >> A >> B;
        cout << floor_sum(N, A, B, M) << endl;
    }
}

int main() {
    //mini_test();
    YosupoSumOfFloorOfLinear();
}
