#define PROBLEM "https://judge.yosupo.jp/problem/aplusb"
#define ERROR 1e-8

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "simplex_method.hpp"

void my_test() {
    // A : n×m行列, n <= mでないと止まらない
    // min cx s.t. Ax = b, x >= 0
    DoubleMatrix<double> A(5, 16, 0);
    vector<double> b(5, 0), c(16, 0);
    for (int i = 0; i < 16; ++i) A[0][i] = 1;
    b = {100, 70, 75, 80, 85};
    for (int i = 0; i < 4; ++i) {
        for (int bit = 0; bit < 16; ++bit) {
            if (bit & (1 << i)) A[i+1][bit] = 1;
        }
    }
    for (int bit = 0; bit < 16; ++bit) {
        if (__builtin_popcount(bit) == 1) c[bit] = -1;
    }
    
    vector<double> res;
    STATE state;
    double opt = -Simplex(A, b, c, res, state);
    
    // 70, 75, 80, 85
    // (2)：5
    // (3)：10
    // (4)：15
    // (1, 2, 3, 4)：70
    double correct_opt = 30;
    vector<double> correct_res = {0, 0, 5, 0, 10, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0, 70};
    assert(abs(opt - correct_opt) <= ERROR);
    for (int i = 0; i < res.size(); ++i) assert(abs(res[i] - correct_res[i]) <= ERROR);
}

// SRM 231 Mixture
class Mixture {
public:
    long double mix(vector<int> mix, vector<string> ab) {
        vector<vector<int>> vab(ab.size());
        for (int i = 0; i < ab.size(); ++i) {
            istringstream sin(ab[i]);
            vector<int> tmp;
            int num;
            while (sin >> num) tmp.push_back(num);
            vab[i] = tmp;
        }
        int n = vab[0].size()-1, m = vab.size();
        DoubleMatrix<long double> A(n, m);
        vector<long double> b(n), c(m);
        
        for (int i = 0; i < n; ++i) b[i] = mix[i];
        for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) A[i][j] = vab[j][i];
        for (int i = 0; i < m; ++i) c[i] = vab[i][n];
        
        vector<long double> res;
        STATE state;
        long double opt = Simplex(A, b, c, res, state);
        
        if (state == OPTIMIZED) return opt;
        else return -1;
    }
};

void librarychecker() {
    int A, B;
    cin >> A >> B;
    cout << A + B << endl;
}

int main() {
    my_test();
    librarychecker();
}
