//
// シェルソート
//
// cf.
//   ソートを極める！ 〜 なぜソートを学ぶのか 〜
//     https://qiita.com/drken/items/44c60118ab3703f7727f
//

#include <iostream>
#include <vector>
using namespace std;

const int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1};

void ShellSort(vector<int> &a) {
    for (int i = 0; i < 8; ++i) {
        int g = gaps[i];

        /* ギャップ g で挿入ソート */
        for (int i = g; i < (int)a.size(); ++i) {
            int v = a[i]; // 挿入したい値
            int j = i;    // 挿入する場所
            for (; j >= g; j -= g) {
                if (a[j - g] > v) {
                    a[j] = a[j - g];
                }
                else break;
            }
            a[j] = v;
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

    /* シェルソート */
    ShellSort(a);

    return 0;
}
