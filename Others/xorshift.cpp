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
