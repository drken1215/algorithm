//
// ランレングス圧縮
//
// verified:
//   第17回 アルゴリズム実技検定(過去問) E - 連長圧縮
//     https://atcoder.jp/contests/past17-open/tasks/past17_e
//


#include <bits/stdc++.h>
using namespace std;


// ランレングス圧縮
template<class Str = string> vector<pair<char,int>> run_length(const Str &S) {
    vector<pair<char,int>> res;
    for (int i = 0; i < (int)S.size();) {
        int j = i;
        while (j < (int)S.size() && S[j] == S[i]) j++;
        res.emplace_back(S[i], j - i);
        i = j;
    }
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// 第17回 アルゴリズム実技検定(過去問) E - 連長圧縮
void PAST_17_E() {
    string S;
    cin >> S;
    auto res = run_length(S);
    for (auto [c, num] : res) cout << c << " " << num << " ";
    cout << endl;
}


int main() {
    PAST_17_E();
}