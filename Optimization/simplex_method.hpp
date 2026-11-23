#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "double_matrix.hpp"

// 2 段階単体法
template<class T> int PreProcess(DoubleMatrix<T> &A, vector<T> &b) {
    const long double INF = 1LL<<60;
    const long double EPS = 1e-10;
    int rank = 0;
    DoubleMatrix<T> M(A.size(), A[0].size()+1);
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
template<class T> T Simplex(DoubleMatrix<T> A, vector<T> b, vector<T> c, vector<T> &res, STATE &state) {
    const long double INF = 1LL<<60;
    const long double EPS = 1e-10;
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
            A.resize(n+1);
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
