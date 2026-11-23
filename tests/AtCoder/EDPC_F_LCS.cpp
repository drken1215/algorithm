#define PROBLEM "https://atcoder.jp/contests/dp/tasks/dp_f"

#include <bits/stdc++.h>
using namespace std;

#include "lcs.hpp"

int main() {
    string S, T;
    cin >> S >> T;
    cout << find_lcs(S, T) << endl;
}
