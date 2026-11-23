#define PROBLEM "https://judge.yosupo.jp/problem/aplusb"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "xorshift.hpp"

void my_test() {
    // N 回サイコロを振って揺らぎを見る
    const int N = 100000;
    const double stdev = sqrt(N * 5) / 6;
    xorshift rng;
    vector<int> con(6, 0);
    for (int i = 0; i < N; ++i) {
        con[rng.rand_int(0, 6)]++;
    }
    for (int i = 0; i < 6; ++i) {
        assert(abs(con[i] - N / 6) <= stdev * 5);
    }
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
