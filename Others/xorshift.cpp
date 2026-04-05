//
// XorShift
//
// cf.
//   ビット演算 (bit 演算) の使い方を総特集！ 〜 マスクビットから bit DP まで 〜
//     https://qiita.com/drken/items/7c6ff2aa4d8fce1c9361
//


#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;


// xor128による乱数生成、周期は2^128-1
unsigned int rand_int() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}
int rand_int(int minv, int maxv) {
    return rand_int() % (maxv - minv + 1) + minv;
}
long long rand_ll(long long minv, long long maxv) {
    long long a = rand_int(), b = rand_int();
    return (a * (1LL<<29) + b) % (maxv - minv + 1) + minv;
}
template<class T> void shuffle(vector<T>& vec) {
    int n = vec.size();
    for (int i = n - 1; i > 0; --i) {
        int k = rand_int() % (i + 1);
        swap(vec[i], vec[k]);
    }
}


//------------------------------//
// Examples
//------------------------------//

int main() {
    // 100,000 回サイコロを振って、それぞれの目が出る回数を数えてみる
    int count[6] = {0};
    for (int i = 0; i < 100000; ++i) {
        int rng = rand_int() % 6;
        count[rng]++;
    }
    for (int i = 0; i < 6; ++i) {
        cout << i+1 << ": " << count[i] << " 回" << endl;
    }
    return 0;
}
