#include <bits/stdc++.h>
using namespace std;

// G[v] := 頂点 v の行き先の頂点
vector<int> detect_cycle(const vector<int> &G) {
    // 頂点を 1 つ好きに選んで、N 回移動する
    int v = 0;
    for (int i = 0; i < (int)G.size(); ++i) v = G[v];
    
    // 得られた頂点から出発してサイクルを検出する
    vector<int> res;
    int start = v;
    do {
        res.push_back(v);
        v = G[v];
    } while (v != start);
    return res;
}

int main() {
    int N;
    cin >> N;
    vector<int> A(N);
    for (int i = 0; i < N; ++i) { cin >> A[i]; --A[i]; }
    
    // 閉路検出
    const auto &res = detect_cycle(A);
    cout << res.size() << endl;
    for (auto v : res) cout << v+1 << " ";
    cout << endl;
}
