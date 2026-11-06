//
// ボゴソート
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
#include <cstdlib>
using namespace std;

void Shuffle(vector<int> &a) {
    int n = (int)a.size();
    for (int i = n - 1; i > 0; --i) {
        int k = rand() % (i + 1);
        swap(a[i], a[k]);
    }
}

void BogoSort(vector<int> &a) {
    while (true) {
        bool swapped = true;
        for (int i = 0; i + 1 < (int)a.size(); ++i) {
            if (a[i] > a[i + 1]) {
                swapped = false;
            }
        }
        if (swapped) break;
        Shuffle(a);
    }
}

int main() {
    int n; // 要素数
    cin >> n;
    vector<int> a(n); // 整列したい配列ベクトル (サイズ を n に初期化)
    for (int i = 0; i < n; ++i) {
        cin >> a[i]; // 整列したい配列を取得
    }

    /* ボゴソート */
    BogoSort(a);

    return 0;
}
