//
// Merge Sort
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

/* 配列 a の [left, right) をソートします */
void MergeSort(vector<int> &a, int left, int right) {
    if (right - left == 1) return;
    int mid = left + (right - left) / 2;

    // 左半分 [left, mid) をソート
    MergeSort(a, left, mid);

    // 右半分 [mid, right) をソート
    MergeSort(a, mid, right);

    // 一旦「左」と「右」のソート結果をコピーしておく (右側は左右反転)
    vector<int> buf;
    for (int i = left; i < mid; ++i) buf.push_back(a[i]);
    for (int i = right-1; i >= mid; --i) buf.push_back(a[i]);

    // マージする
    int iterator_left = 0;                    // 左側のイテレータ
    int iterator_right = (int)buf.size() - 1; // 右側のイテレータ
    for (int i = left; i < right; ++i) {
        // 左側採用
        if (buf[iterator_left] <= buf[iterator_right]) {
            a[i] = buf[iterator_left++];
        }
        // 右側採用
        else {
            a[i] = buf[iterator_right--];
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

    /* マージソート */
    MergeSort(a, 0, n);

    return 0;
}
