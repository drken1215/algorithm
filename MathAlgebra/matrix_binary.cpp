//
// binary 行列 (行列累乗と、掃き出し法)
//
// verified:
//   みんなのプロコン 2019 E - Odd Subrectangles
//     https://atcoder.jp/contests/yahoo-procon2019-qual/tasks/yahoo_procon2019_qual_e
// 


#include <iostream>
#include <vector>
#include <bitset>
using namespace std;


const int MAX_ROW = 510; // to be set
const int MAX_COL = 510; // to be set
struct BitMatrix {
    int n, m;
    bitset<MAX_COL> val[MAX_ROW];
    BitMatrix(int n_ = 1, int m_ = 1) {n = n_; m = m_;}
    inline bitset<MAX_COL>& operator [] (int i) {return val[i];}
    inline friend ostream& operator << (ostream& s, BitMatrix M) {
        s << endl; 
        for (int i = 0; i < M.n; ++i) {
            for (int j = 0; j < M.m; ++j) s << M.val[i][j];
            s << endl;
        }
        return s;
    }
};

inline BitMatrix operator * (BitMatrix A, BitMatrix B) {
    BitMatrix R(A.n, B.m);
    BitMatrix tB(B.m, B.n);
    for (int i = 0; i < tB.n; ++i) for (int j = 0; j < tB.m; ++j) tB[i][j] = B[j][i];
    for (int i = 0; i < R.n; ++i) for (int j = 0; j < R.m; ++j) R[i][j] = (A[i] & tB[j]).any();
    return R;
}

inline BitMatrix pow(BitMatrix A, long long n) {
    BitMatrix R(A.n, A.n);
    for (int i = 0; i < A.n; ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * A;
        A = A * A;
        n >>= 1;
    }
    return R;
}

int calc_rank(BitMatrix &A) {
    int r = 0;
	for (int i = 0; i < A.m; ++i) {
		int pivot = -1;
		for (int j = r; j < A.n; ++j) {
			if (A[j][i]) {
				pivot = j;
				break;
			}
		}
		if (pivot != -1) {
			swap(A[pivot], A[r]);
			for (int j = 0; j < A.n; ++j) if (j != r && A[j][i]) A[j] ^= A[r];
			++r;
		}
	}
    return r;
}

vector<vector<int> > linear_equation(BitMatrix A, vector<int> b) {
	int rank = 0;
    for (int i = 0; i < A.n; ++i) { A[i][A.m] = b[i]; }
    
    vector<int> core, rem;
	for (int i = 0; i < A.m; ++i) {
		int pivot = -1;
		for (int j = rank; j < A.n; ++j) {
			if (A[j][i]) {
				pivot = j;
				break;
			}
		}
		if (pivot != -1) {
            core.push_back(i);
			swap(A[pivot], A[rank]);
			for (int j = 0; j < A.n; ++j) if (j != rank && A[j][i]) A[j] ^= A[rank];
			++rank;
		}
        else rem.push_back(i);
	}
    
    vector<vector<int> > res;
    for (int i = rank; i < A.n; ++i) 
        if (A[i][A.m]) return res;     // return -1;
    
    vector<int> sol(A.m, 0);
    for (int i = 0; i < core.size(); ++i) sol[core[i]] = A[i][A.m];
    res.push_back(sol);
    
    for (int i = 0; i < rem.size(); ++i) {
        vector<int> temp(A.m, 0);
        temp[rem[i]] = 1;
        for (int j = 0; j < core.size(); ++j) temp[core[j]] = A[j][rem[i]];
        res.push_back(temp);
    }
    
    return res;     // return A[0].size()-rank;
};


const int MOD = 998244353;
long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

int main() { 
    int N, M; cin >> N >> M;
    BitMatrix A(N, M);
    for (int i = 0; i < N; ++i) for (int j = 0; j < M; ++j) {
            int a; cin >> a;
            if (a) A[i][j] = 1;
        }
    int r = calc_rank(A);  
    cout << (modpow(2LL, N+M-1, MOD) - modpow(2LL, N+M-r-1, MOD) + MOD) % MOD << endl;
}
