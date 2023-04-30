//
// 計数ソート
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

const int MAX = 100000; // 配列の値は 100000 未満だとします

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* 各要素の個数をカウントします */
    /* num[v] := v の個数 */
    int num[MAX] = {0};
    for (int i = 0; i < n; ++i) {
        ++num[a[i]]; // a[i] をカウントします
    }

    /* num の累積和をとります */
    /* sum[v] := v 以下の値の個数 */
    int sum[MAX] = {0};
    for (int v = 1; v < MAX; ++v) {
        sum[v] = sum[v-1] + num[v];
    }

    /* sum を元にソート */
    /* sorted: a をソートしたもの */
    vector<int> sorted(n);
    for (int i = n-1; i >= 0; --i) {
        sorted[--sum[a[i]]] = a[i];
    }

    return 0;
}
