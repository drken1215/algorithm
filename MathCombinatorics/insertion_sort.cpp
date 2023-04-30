//
// Insertion Sort
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* 挿入前の配列を出力してみる */
    cout << "Before: ";
    for (int i = 0; i < n; ++i) cout << a[i] << " ";
    cout << endl;

    /* 挿入ソート */
    for (int i = 1; i < n; ++i) {
        int v = a[i]; // 挿入したい値

        // v を挿入する適切な場所 j を探す
        int j = i;
        for (; j > 0; --j) {
            if (a[j-1] > v) { // v より大きいやつは 1 個後ろに移す
                a[j] = a[j-1];
            }
            else break; // v 以下になったら止める
        }
        a[j] = v; // 最後に j 番目に v を挿入

        /* 各ステップの配列を出力してみる */
        cout << "After Step " << i << ": ";
        for (int i = 0; i < n; ++i) cout << a[i] << " ";
        cout << endl;
    }

    return 0;
}
