//
// Zobrist Hash (集合に対するハッシュ)
//
// reference:
//   https://codeforces.com/blog/entry/62393
//
// verified:
//   AtCoder ABC 250 E - Prefix Equality
//     https://atcoder.jp/contests/abc250/tasks/abc250_e
//


#include <bits/stdc++.h>
using namespace std;


// 整数値 x にハッシュ値を割り当てる関数
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    uint64_t operator() (uint64_t x) const {
        static const uint64_t FIXED_RANDOM =
            chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
} rng;


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ABC_250_E() {
    // Zobrist Hash
    // hashA[i] := 数列 a の先頭から i 個分の集合を表すハッシュ値
    int N;
    cin >> N;
    
    // N 個の入力を受け取って、Zobrist Hash を返す関数
    auto hashing = [&]() -> vector<uint64_t> {
        vector<uint64_t> hash(N + 1, 0);
        set<long long> S;  // すでにあるかを確認する
        
        for (int i = 0; i < N; ++i) {
            long long x;
            cin >> x;
            
            // すでに含まれている場合は何もしない
            if (S.count(x)) {
                hash[i + 1] = hash[i];
                continue;
            }
            S.insert(x);
            hash[i + 1] = hash[i] ^ rng(x);
        }
        return hash;
    };
    
    const auto& hashA = hashing();
    const auto& hashB = hashing();
    
    // クエリ処理
    int Q;
    cin >> Q;
    while (Q--) {
        int x, y;
        cin >> x >> y;
        if (hashA[x] == hashB[y]) cout << "Yes" << endl;
        else cout << "No" << endl;
    }
}

int main() {
    ABC_250_E();
}


