//
// コムソート
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

void CombSort(vector<int> &a) {
    int g = (int)a.size();
    while (true) {
        bool nonswaped = true;      // swap されてないかどうか
        g = max((g * 10) / 13, 1);  // gap を縮める (1 以上にはする)

        /* 間隔 g でバブルソートっぽく */
        for (int i = 0; i + g < (int)a.size(); ++i) {
            if (a[i] > a[i + g]) {
                swap(a[i], a[i + g]);
                nonswaped = false;
            }
        }
        if (g == 1 && nonswaped) break; // g = 1 で swap なしなら break
    }
}

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* コムソート */
    CombSort(a);

    return 0;
}
