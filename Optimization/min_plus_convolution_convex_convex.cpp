//
// Min-Plus Convolution (Convex and Convex), O(N + M)
//   区分線形凸関数の傾き (階差数列) の列をマージしてソートすればよい
//
// verified:
//   Library Checker - Min Plus Convolution (Convex and Convex)
//     https://judge.yosupo.jp/problem/min_plus_convolution_convex_convex
//
// Reference:
//   AtCoder Lecture Note: 凸数列の min-plus 畳み込み
//     https://info.atcoder.jp/entry/algorithm_lectures/min_plus_convolution
//   


#include <bits/stdc++.h>
using namespace std;


// Min-Plus Convolution (Convex and Convex)
template<class T> vector<T> MinPlusConvolutionConvexConvex(vector<T> &A, vector<T> &B) {
    int N = (int)A.size(), M = (int)B.size();
    if (N == 0 || M == 0) return {};
    vector<T> res(N + M - 1);
    int a = 0, b = 0;
    res[0] = A[0] + B[0];
    for (int i = 1; i < N + M - 1; i++) {
        if (b == M - 1 || (a < N - 1 && A[a + 1] + B[b] < A[a] + B[b + 1])) {
            res[i] = A[++a] + B[b];
        } else {
            res[i] = A[a] + B[++b];
        }
    }
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Min Plus Convolution (Convex and Convex)
void LibraryCheckerMinPlusConvolutionConvexConvex() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    long long N, M;
    cin >> N >> M;
    vector<long long> A(N), B(M);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < M; i++) cin >> B[i];
    auto res = MinPlusConvolutionConvexConvex(A, B);
    for (int i = 0; i < (int)res.size(); i++) cout << res[i] << " ";
    cout << endl;
}


int main() {
    LibraryCheckerMinPlusConvolutionConvexConvex();
}