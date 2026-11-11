#define PROBLEM "https://yukicoder.me/problems/no/704"

#include <bits/stdc++.h>
using namespace std;

#include "monge_shortest_path_by_cht.hpp"

int main() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx + dy;
    };
    auto res = MongeShortestPathByCHT<long long>(N, func);
    cout << res[N].first << endl;
}
