//
// nCk mod m、単純計算
//
// cf.
//   よくやる二項係数 (nCk mod. p)、逆元 (a^-1 mod. p) の求め方
//     http://drken1215.hatenablog.com/entry/2018/06/08/210000
//
// verified:
//   CODECHEF July Challenge 2018 Division 2 D - No Minimum No Maximum
//     https://www.codechef.com/problems/NMNMX
//


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


inline long long mod(long long a, long long m) {
    return (a % m + m) % m;
}

inline long long inv(long long a, long long m) {
    long long b = m, u = 1, v = 0;
    while (b) {
        long long t = a/b;
        a -= t*b; swap(a, b);
        u -= t*v; swap(u, v);
    }
    return mod(u, m);
}

const int MAX_R = 510000;
long long DCom[MAX_R];

void calc_Dcom(long long n, long long p) {
    long long t = 1; DCom[0] = t;
    for (int i = 1; i < MAX_R; ++i) {
        t = mod(t * (n-i+1), p);
        t = mod(t * inv(i, p), p);
        DCom[i] = t;
    }
}



//------------------------------//
// Examples
//------------------------------//

int main() {

}
