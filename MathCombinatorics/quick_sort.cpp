//
// Quick Sort
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

/* 配列 a の [left, right) をソートします */
void QuickSort(vector<int> &a, int left, int right) {
    if (right - left <= 1) return;

    int pivot_index = (left + right) / 2;  // 適当にここでは中点とします
    int pivot = a[pivot_index];
    swap(a[pivot_index], a[right - 1]);    // pivot と右端を swap

    int i = left; // iterator
    for (int j = left; j < right - 1; ++j) { // j は全体を眺めて
        if (a[j] < pivot) { // pivot 未満のがあったら左に詰めていく
            swap(a[i++], a[j]);
        }
    }
    swap(a[i], a[right - 1]); // pivot を適切な場所に挿入

    /* 再帰的に解く */
    QuickSort(a, left, i);    // 左半分 (pivot 未満)
    QuickSort(a, i + 1, right); // 右半分 (pivot 以上)
}

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* クイックソート */
    QuickSort(a, 0, n);

    return 0;
}
