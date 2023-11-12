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

const long double INF = 1LL<<60;
const long double EPS = 1e-10;

template<class T> int PreProcess(Matrix<T> &A, vector<T> &b) {
    int rank = 0;
    Matrix<T> M(A.size(), A[0].size()+1);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) M[i][j] = A[i][j];
        M[i][A[0].size()] = b[i];
    }
    
    for (int i = 0; i < A[0].size(); ++i) {
        int pivot = rank;
        T max = 0;
        for (int j = rank; j < M.size(); ++j) {
            if (max < abs(M[j][i])) {
                pivot = j;
                max = abs(M[j][i]);
            }
        }
        if (max > EPS) {
            swap(M[pivot], M[rank]);
            T div = M[rank][i];
            for (int j = 0; j < M[0].size(); ++j) {
                M[rank][j] = M[rank][j] / div;
            }
            for (int j = 0; j < M.size(); ++j) {
                if (j != rank && abs(M[j][i]) > EPS) {
                    T fac = M[j][i];
                    for (int k = 0; k < M[0].size(); ++k) {
                        M[j][k] = M[j][k] - M[rank][k] * fac;
                    }
                }
            }
            ++rank;
        }
    }
    A.resize(rank);
    b.resize(rank);
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < A[0].size(); ++j) A[i][j] = M[i][j];
        b[i] = M[i].back();
    }
    return rank;
}

enum STATE { OPTIMIZED, INFEASIBLE, UNBOUNDED };
template<class T> T Simplex(Matrix<T> A, vector<T> b, vector<T> c, vector<T> &res, STATE &state) {
    res = vector<T>();
    
    // let A be row-fullrank
    PreProcess(A, b);
    
    // build tableau
    int n = A.size(), m = A[0].size();
    for (int i = 0; i < n; ++i) if (b[i] < 0) {
        b[i] *= -1;
        for (int j = 0; j < m; ++j) A[i][j] *= -1;
    }
    vector<int> base(n), non(m);
    for (int i = 0; i < n; ++i) base[i] = m+i;
    for (int j = 0; j < m; ++j) non[j] = j;
    
    A.resize(n+2, n+m+1);
    for (int i = 0; i < n; ++i) A[i][m+i] = 1, A[i][n+m] = b[i];
    for (int i = 0; i < n; ++i) {
        A[n][n+m] += A[i][n+m];
        for (int j = 0; j < m; ++j) A[n][j] += A[i][j];
    }
    for (int j = 0; j < m; ++j) A[n+1][j] = -c[j];
    
    // start simplex
    for (int phase = 0; phase < 2; ++phase) {
        while (true) {
            int nn = -1, nb = -1;
            for (int i = 0; i < non.size(); ++i) {
                // We cannot let slack value move to the base
                if (non[i] >= m) continue;
                if (A[n][non[i]] > EPS) {
                    if (nn == -1) nn = i;
                    else if (non[i] < non[nn]) {
                        // Bland's smallest subscript rule
                        nn = i;
                    }
                }
            }
            if (nn == -1) {
                if (phase == 1) {
                    // All done!
                    break;
                }
                if (A[n][A[0].size()-1] > EPS) {
                    // No feasible solution!
                    state = INFEASIBLE;
                    return -1;
                }
                
                // detail postprocess of phase 0
                bool ok = true;
                for (int i = 0; i < base.size(); ++i) {
                    if (base[i] >= m) {
                        // If base doesn't contain slack, go to phase 1
                        ok = false;
                        nb = i;
                        break;
                    }
                }
                if (ok) break;
                for (int i = 0; i < non.size(); ++i) {
                    if (non[i] < m && abs(A[nb][non[i]]) > EPS) {
                        // slack has base, continue phase 0
                        nn = i;
                        break;
                    }
                }
                if (nn == -1) {
                    // It can't be happened (A is row-fullrank)
                    break;
                }
            }
            int col = non[nn];
            
            if (nb == -1) {
                T min_ratio = INF;
                for (int i = 0; i < base.size(); ++i) if (A[i][col] > EPS) {
                    T tmp = A[i][A[0].size()-1] / A[i][col];
                    if (min_ratio > tmp + EPS) {
                        min_ratio = tmp, nb = i;
                    } else if (min_ratio >= tmp - EPS) {
                        if (nb == -1) nb = i;
                        else if (base[i] < base[nb]) {
                            // Bland's smallest subscript rule
                            nb = i;
                        }
                    }
                }
                if (nb == -1) {
                    // It cannot happen at the 1st stage
                    state = UNBOUNDED;
                    return -1;
                }
            }
            int row = nb;
            
            T piv = A[row][col];
            for (int j = 0; j < A[0].size(); ++j) A[row][j] = A[row][j] / piv;
            for (int i = 0; i < A.size(); ++i) {
                if (i == row) continue;
                T pn = A[i][col];
                for (int j = 0; j < A[0].size(); ++j) A[i][j] = A[i][j] - A[row][j] * pn;
            }
            swap(base[nb], non[nn]);
        }
        
        if (phase == 0) {
            swap(A[n], A[n+1]);
            A.val.pop_back();
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < A.size(); ++j)
                    A[j].erase(A[j].begin() + m);
            
            for (int i = 0; i < non.size(); ++i)
                if (non[i] >= m)
                    non.erase(non.begin() + i--);
        }
    }
    res = vector<T>(m, 0);
    for (int i = 0; i < base.size(); ++i) res[base[i]] = A[i].back();
    state = OPTIMIZED;
    return A[n].back();
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void small_test() {
    // A : n×m行列, n <= mでないと止まらない
    // min cx s.t. Ax = b, x >= 0
    Matrix<double> A(5, 16, 0);
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
    double opt = Simplex(A, b, c, res, state);
    
    // 70, 75, 80, 85
    // (2)：5
    // (3)：10
    // (4)：15
    // (1, 2, 3, 4)：70
    cout << -opt << endl;
    for (int i = 0; i < res.size(); ++i) cout << res[i] << " ";
    cout << endl;
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
        Matrix<long double> A(n, m);
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


int main() {
    small_test();
}
