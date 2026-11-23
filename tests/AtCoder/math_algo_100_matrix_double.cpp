#define PROBLEM "https://atcoder.jp/contests/math-and-algorithm/tasks/math_and_algorithm_bv"
#define ERROR 1e-8

#include <bits/stdc++.h>
using namespace std;

#include "double_matrix.hpp"

int main() {
    using DD = long double;
    int Q;
    cin >> Q;
    while (Q--) {
        DD X, Y, Z;
        long long T;
        cin >> X >> Y >> Z >> T;

        DoubleMatrix<DD> A(3, 3, 0);
        A[0][0] = 1.0 - X, A[0][1] = Y;
        A[1][1] = 1.0 - Y, A[1][2] = Z;
        A[2][2] = 1.0 - Z, A[2][0] = X;
        auto P = A.pow(T);
        vector<DD> ini({1, 1, 1});
        auto res = P * ini;
        cout << fixed << setprecision(10);
        cout << res[0] << " " << res[1] << " " << res[2] << endl;
    }
}
