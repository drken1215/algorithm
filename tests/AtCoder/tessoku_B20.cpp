#define PROBLEM "https://atcoder.jp/contests/tessoku-book/tasks/tessoku_book_cs"

#include <bits/stdc++.h>
using namespace std;

#include "edit_distance.hpp"

int main() {
    string S, T;
    cin >> S >> T;
    cout << find_edit_distance(S, T) << endl;
}
