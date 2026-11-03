//
// 与えられた部分集合の部分集合を列挙
//
// cf.
//   ビット演算 (bit 演算) の使い方を総特集！ 〜 マスクビットから bit DP まで 〜
//     https://qiita.com/drken/items/7c6ff2aa4d8fce1c9361
//


#include <iostream>
#include <vector>
using namespace std;


int main() {
    int A = (1<<2) | (1<<3) | (1<<5) | (1<<7); // {A = {2, 3, 5, 7}
    
    for (int bit = A; ; bit = (bit-1) & A) {
        /* ここに処理を書く */
        
        
        // 最後の 0 で break
        if (!bit) break;
    }
}
