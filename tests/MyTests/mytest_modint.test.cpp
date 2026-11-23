#define PROBLEM "https://judge.yosupo.jp/problem/aplusb"

#include <bits/stdc++.h>
using namespace std;

#include "modint.hpp"

void my_test() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    auto check = [&](mint val, mint res) -> void {
        //cout << val << endl;
        assert(val == res);
    };
    
    // arithmetic operation
    mint a = -3;
    mint b = 5;
    check(a + b, 2);
    check(a - b, 998244345);  // -8
    check(a * b, 998244338);  // -15
    a += b;
    check(a, 2);
    a -= 3;
    check(a, 998244352);  // -1
    a *= 6;
    check(a, 998244347);  // -6
    a /= 2;
    check(a, 998244350);  // -3

    // increment
    check(++a, 998244351);  // -2
    check(a, 998244351);  // -2
    check(--a, 998244350);  // -3
    check(a, 998244350);  // -3
    check(b++, 5);
    check(b, 6);
    check(b--, 6);
    check(b, 5);
    
    // division
    a = 6, b = 2;
    check(a / b, 3);
    mint c = b / a;
    check(c, 332748118);
    check(c * a, 2);
    
    // pow, inv
    check(b.pow(20), 1048576);
    check(b.inv(), 499122177);
}

void librarychecker() {
    int A, B;
    cin >> A >> B;
    cout << A + B << endl;
}

int main() {
    my_test();
    librarychecker();
}
