#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/1328"

#include <bits/stdc++.h>
using namespace std;

#include "double_matrix.hpp"

int main() {
    using D = long double;
    
    auto dpow = [&](D a, int n) -> D {
        D res = 1.0;
        for (int i = 0; i < n; ++i) res *= a;
        return res;
    };
    auto func = [&](const vector<D> &coef, int i) -> D {
        D res = 0.0;
        for (int p = 0; p < (int)coef.size(); ++p)
            res += coef[p] * pow(i, p);
        return res;
    };
    
    int d;
    while (cin >> d, d) {
        vector<D> v(d + 3);
        for (int i = 0; i < d + 3; ++i) cin >> v[i];

        bool finish = false;
        int res = 0;
        for (int i = 0; i < d + 3 && !finish; ++i) {
            for (int j = i + 1; j < d + 3 && !finish; ++j) {
                DoubleMatrix<D> A(d + 1, d + 1);
                A.set_eps(1e-5);
                vector<D> b(d + 1);
                for (int k = 0, iter = 0; k < d+3; ++k) {
                    if (k == i || k == j) continue;
                    for (int p = 0; p < d + 1; ++p) {
                        A[iter][p] = pow(k, p);
                        b[iter] = v[k];
                    }
                    ++iter;
                }
                vector<D> ans;
                auto rank = linear_equation(A, b, ans);
                if (ans.empty()) continue;
                D vi = func(ans, i), vj = func(ans, j);
                int num = 0;
                if (fabs(vi - v[i]) > A.EPS) res = i, ++num;
                if (fabs(vj - v[j]) > A.EPS) res = j, ++num;
                if (num == 1) goto end;
            }
        }
        end:
        cout << res << endl;
    }
}
