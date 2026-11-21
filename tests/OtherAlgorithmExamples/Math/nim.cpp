//
// Nim AI (three moutaions)
//   win pattern: xor = 0
//


#include <iostream>
#include <cstdlib>
using namespace std;

int main() {
    // 3 山の石の個数を入力
    int A, B, C;
    cout << "The number of the first pile: ";
    cin >> A;
    cout << "The number of the second pile: ";
    cin >> B;
    cout << "The number of the third pile: ";
    cin >> C;

    // 先手と後手の戦略
    string ans;
    cout << "Are you first? (yes / no): " << endl;
    cin >> ans;
    bool yours = (ans == "yes" ? true : false);

    // 石を取れなくなるまでプレイする
    while (A + B + C > 0) {
        cout << "-----------------" << endl;
        
        // プレイヤーのターン (例外処理は省略)
        if (yours) {
            // 石の山と取る石の個数を入力
            int which, num;
            cout << "Your Turn" << endl;
            cout << "Current state: ("
                 << A << ", "
                 << B << ", "
                 << C << ")" << endl;
            cout << "Which pile? (1 or 2 or 3): ";
            cin >> which;
            cout << "The number of stones? : ";
            cin >> num;

            // 石の個数を変更
            if (which == 1) A -= num;
            else if (which == 2) B -= num;
            else C -= num;
        }
        // AI のターン
        else {
            cout << "AI's turn" << endl;
            cout << "Current state: ("
                 << A << ", "
                 << B << ", "
                 << C << ")" << endl;

            // A, B, C の XOR 和を S とする
            int S = A ^ B ^ C;

            // S > 0 のときは勝ちパターンにする
            if (S > 0) {
                if ((A ^ S) < A) A ^= S;
                else if ((B ^ S) < B) B ^= S;
                else C ^= S;
            }
            // それ以外の場合はランダムに石を減らす
            else {
                if (A > 0) A = rand() % A;
                else if (B > 0) B = rand() % B;
                else C = rand() % C;
            }

            cout << "After AI's turn: ("
                 << A << ", "
                 << B << ", "
                 << C << ")" << endl;
        }

        // 手番を入れ替え
        yours = !yours;
    }

    if (yours) cout << "You lose..." << endl;
    else cout << "You win!" << endl;
}
