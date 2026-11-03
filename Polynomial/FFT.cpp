//
// FFT (Fast Fourier Transform)
//
// verified:
//   AtCoder Typical Contest 001 C - 高速フーリエ変換
//     https://atc001.contest.atcoder.jp/tasks/fft_c
//


#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


// FFT
namespace FFT {
    using DD = double;
    const DD PI = acosl(-1);

    struct Comp {
        DD real, imag;
        Comp(DD real = 0, DD imag = 0) : real(real), imag(imag) {}
        friend inline ostream& operator << (ostream &s, const Comp &c) {
            return s << '<' << c.real << ',' << c.imag << '>';
        }
        inline Comp operator + (const Comp &c) {
            return {real + c.real, imag + c.imag};
        }
        inline Comp operator - (const Comp &c) {
            return {real - c.real, imag - c.imag};
        }
        inline Comp operator * (const Comp &c) {
            return {real * c.real - imag * c.imag,
                    real * c.imag + imag * c.real};
        }
        inline Comp operator * (DD a) {
            return {real * a, imag * a};
        }
        inline Comp operator / (DD a) {
            return {real / a, imag / a};
        }
    };

    // FFT
    void trans(vector<Comp> &v, bool inv = false) {
        int n = (int)v.size();
        for (int i = 0, j = 1; j < n-1; j++) {
            for (int k = n>>1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int t = 2; t <= n; t <<= 1) {
            DD ang = acosl(-1.0) * 2 / t;
            if (inv) ang = -ang;
            for (int i = 0; i < n; i += t) {
                for (int j = 0; j < t/2; ++j) {
                    Comp w = {cos(ang * j), sin(ang * j)};
                    int j1 = i + j, j2 = i + j + t/2;
                    Comp c1 = v[j1], c2 = v[j2] * w;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2;
                }
            }
        }
        if (inv) for (int i = 0; i < n; ++i) v[i] = v[i]/n;
    }
    
    // A * B
    vector<long long> mult(const vector<long long> &A,
                           const vector<long long> &B) {
        int size_a = 1; while (size_a < A.size()) size_a <<= 1;
        int size_b = 1; while (size_b < B.size()) size_b <<= 1;
        int size_fft = max(size_a, size_b) << 1;
        
        vector<Comp> cA(size_fft), cB(size_fft), cC(size_fft);
        for (int i = 0; i < A.size(); ++i) cA[i] = {(DD)A[i], 0};
        for (int i = 0; i < B.size(); ++i) cB[i] = {(DD)B[i], 0};
        
        trans(cA); trans(cB);
        for (int i = 0; i < size_fft; ++i) cC[i] = cA[i] * cB[i];
        trans(cC, true);
        
        vector<long long> res((int)A.size() + (int)B.size() - 1);
        for (int i = 0; i < res.size(); ++i) {
            res[i] = (long long)(cC[i].real + 0.5);
        }
        return res;
    }
};



//------------------------------//
// Examples
//------------------------------//

int main() {
    int N;
    while (cin >> N) {
        vector<long long> a(N), b(N);
        for (int i = 0; i < N; ++i) cin >> a[i] >> b[i];
        auto res = FFT::mult(a, b);
        cout << 0 << endl;
        for (int i = 0; i < N*2-1; ++i) cout << res[i] << endl;
    }
}
