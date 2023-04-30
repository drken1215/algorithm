//
// Selection Sort
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

void SelectionSort(vector<int> &a) {
    for (int i = 0; i < (int)a.size(); ++i) {
        int min_index = i; // i 番目以降で一番小さい値を探す
        for (int j = min_index; j < (int)a.size(); ++j) {
            if (a[j] < a[min_index]) {
                min_index = j; // より小さいのがあったら更新
            }
        }
        // 交換
        swap(a[i], a[min_index]);
    }
}

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* 選択ソート */
    SelectionSort(a);

    return 0;
}
