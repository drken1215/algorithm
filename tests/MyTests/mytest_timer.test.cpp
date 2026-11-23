#define PROBLEM "https://judge.yosupo.jp/problem/aplusb"

#include <bits/stdc++.h>
using namespace std;

#include "timer.hpp"

void my_test() {
    Timer timer1, timer2(true);
    timer1.start();
    int iter = 0;
    while (timer1() < 0.001) iter++;
    double t1 = timer1(), t2 = timer2();
    assert(abs(t1 - t2) < 0.1);
    timer1.stop();
    while (timer2() < 0.02) iter++;
    timer1.start();
    t1 = timer1(), t2 = timer2();
    assert(t2 > t1 + 0.019);
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
