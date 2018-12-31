//
// FFT (Fast Fourier Transform)
//
// cf.
//
//
// verified:
//   AtCoder Typical Contest 001 C - 高速フーリエ変換
//     https://atc001.contest.atcoder.jp/tasks/fft_c
//

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;


// complex number
struct ComplexNumber {
    double real, imag;
    inline ComplexNumber& operator = (const ComplexNumber &c) {real = c.real; imag = c.imag; return *this;}
    friend inline ostream& operator << (ostream &s, const ComplexNumber &c) {return s<<'<'<<c.real<<','<<c.imag<<'>';}
};
inline ComplexNumber operator + (const ComplexNumber &x, const ComplexNumber &y) {
    return {x.real + y.real, x.imag + y.imag};
}
inline ComplexNumber operator - (const ComplexNumber &x, const ComplexNumber &y) {
    return {x.real - y.real, x.imag - y.imag};
}
inline ComplexNumber operator * (const ComplexNumber &x, const ComplexNumber &y) {
    return {x.real * y.real - x.imag * y.imag, x.real * y.imag + x.imag * y.real};
}
inline ComplexNumber operator * (const ComplexNumber &x, double a) {
    return {x.real * a, x.imag * a};
}
inline ComplexNumber operator / (const ComplexNumber &x, double a) {
    return {x.real / a, x.imag / a};
}

// fft
namespace FFT {
    void trans(vector<ComplexNumber> &v, bool inv = false) {
        int n = (int)v.size();
        for (int i = 0, j = 1; j < n-1; j++) {
            for (int k = n>>1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int t = 2; t <= n; t <<= 1) {
            double ang = acos(-1.0) * 2 / t;
            if (inv) ang = -ang;
            for (int i = 0; i < n; i += t) {
                for (int j = 0; j < t/2; ++j) {
                    ComplexNumber w = {cos(ang * j), sin(ang * j)};
                    int j1 = i + j, j2 = i + j + t/2;
                    ComplexNumber c1 = v[j1], c2 = v[j2] * w;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2;
                }
            }
        }
        if (inv) for (int i = 0; i < n; ++i) v[i] = v[i]/n;
    }
    
    // C is A*B
    vector<long long> mult(vector<long long> A, vector<long long> B) {
        int size_a = 1; while (size_a < A.size()) size_a <<= 1;
        int size_b = 1; while (size_b < B.size()) size_b <<= 1;
        int size_fft = max(size_a, size_b) << 1;
        
        vector<ComplexNumber> cA(size_fft), cB(size_fft), cC(size_fft);
        for (int i = 0; i < A.size(); ++i) cA[i] = {(double)A[i], 0};
        for (int i = 0; i < B.size(); ++i) cB[i] = {(double)B[i], 0};
        
        trans(cA); trans(cB);
        for (int i = 0; i < size_fft; ++i) cC[i] = cA[i] * cB[i];
        trans(cC, true);
        
        vector<long long> res((int)A.size() + (int)B.size() - 1);
        for (int i = 0; i < res.size(); ++i) res[i] = (long long)(cC[i].real + 0.5);
        return res;
    }
};



int main() {
    int n; cin >> n;
    vector<long long> A(n), B(n);
    for (int i = 0; i < n; ++i) cin >> A[i] >> B[i];
    
    vector<long long> C = FFT::mult(A, B);
    cout << 0 << endl;
    for (int i = 0; i < n*2-1; ++i) cout << C[i] << endl;
}
