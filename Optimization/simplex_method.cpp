//
// 最小添字規則によるナイーブな2段階単体法
//
// verified:
//   忘れ物の問題 (100 人の生徒のうち、A, B, C, D を忘れた生徒は 70 人, 75 人, 80 人, 85 人)
//  　 1 種類のみを忘れた生徒の人数の最大値を求める
//
//   TopCoder SRM 231 DIV1 Hard Mixture
//     https://community.topcoder.com/stat?c=problem_statement&pm=3942&rd=6520
//


#include <bits/stdc++.h>
using namespace std;


// 単体法（最小添字規則によるナイーブな2段階単体法）, 有理数ライブラリにも適用可能
// A : n×m行列, n <= m でないと止まらない
// min cx s.t. Ax = b, x >= 0
template<class T> struct Matrix {
    vector<vector<T> > val;
    Matrix(int n = 1, int m = 1) {val.clear(); val.resize(n, vector<T>(m));}
    Matrix(int n, int m, T x) {val.clear(); val.resize(n, vector<T>(m, x));}
    void init(int n, int m, T x = 0) {val.clear(); val.resize(n, vector<T>(m, x));}
    void resize(int n) {val.resize(n);}
    void resize(int n, int m, T x = 0) {val.resize(n); for (int i = 0; i < n; ++i) val[i].resize(m, x);}
    int size() {return val.size();}
    inline vector<T>& operator [] (int i) {return val[i];}
    friend ostream& operator << (ostream& s, Matrix<T> M) {s << endl;
        for (int i = 0; i < M.val.size(); ++i) s << M[i] << endl; return s;}
};




//------------------------------//
// Examples
//------------------------------//

void small_test() {

}


/

int main() {
    small_test();
}
