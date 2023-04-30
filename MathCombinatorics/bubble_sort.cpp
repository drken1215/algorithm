//
// Bubble Sort
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

void BubbleSort(vector<int> &a) {
    for (int i = 0; i < (int)a.size() - 1; ++i) {
        for (int j = (int)a.size() - 1; j > i; --j) {
            if (a[j - 1] > a[j]) {  // 大きさが逆転している箇所があったら swap
                swap(a[j - 1], a[j]);
            }
        }
    }
}

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* バブルソート */
    BubbleSort(a);

    return 0;
}
