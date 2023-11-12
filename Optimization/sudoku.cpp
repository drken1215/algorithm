//
// 数独ソルバー
//
// verified:
//   POJ 2676 - Sudoku (solve 9 x 9 sudoku)
//     http://poj.org/problem?id=2676
//
//   AtCoder ABC 327 C - Number Place (judge validity)
//     https://atcoder.jp/contests/abc327/tasks/abc327_c
//


/*
    入力: 以下のいずれかをコンストラクタで与える
    - vector<string> 型変数 (空きマスはデフォルトでは '*')
    - vector<vector<int>> 型変数 (空きマスはデフォルトでは '-1')
 
    // 使用例
    vector<string> input = {
        "53**7****",
        "6**195***",
        "*98****6*",
        "8***6***3",
        "4**8*3**1",
        "7***2***6",
        "*6****28*",
        "***419**5",
        "****8**79"
    };
    Sudoku sudoku(input);
 
    // 数独を解く (すべての解を求める)
    vector<vector<vector<int>>> res = sudoku.solve();
 */


#include <bits/stdc++.h>
using namespace std;


// 未確定を表す数値
const int UNPUT = -1;

// 数独を解くための構造体
class Sudoku {
 public:
    // 盤面を二次元ベクトルで表す
    using Field = vector<vector<int>>;
    
    // コンストラクタ (未確定マスの値を -1 で表す)
    Sudoku(int siz = 3) : siz(siz) {
        init(siz);
    }
    Sudoku(const Sudoku &board)
    : siz(board.siz), field(board.field),
    nums(board.nums), choices(board.choices), is_valid(board.is_valid) {
    }
    Sudoku(const vector<string> &input, char empty_cell = '*') {
        siz = 0;
        while (siz * siz < (int)input.size()) ++siz;
        init(siz);
        assert(input.size() == siz * siz);
        for (int x = 0; x < siz * siz; ++x) {
            assert(input[x].size() == siz * siz);
            for (int y = 0; y < siz * siz; ++y) {
                if (input[x][y] == empty_cell) continue;
                int val = input[x][y] - '0';
                if (!put(x, y, val)) is_valid = false;
            }
        }
    }
    Sudoku(const vector<vector<int>> &input, int empty_cell = -1) {
        siz = 0;
        while (siz * siz < (int)input.size()) ++siz;
        init(siz);
        assert(input.size() == siz * siz);
        for (int x = 0; x < siz * siz; ++x) {
            assert(input[x].size() == siz * siz);
            for (int y = 0; y < siz * siz; ++y) {
                if (input[x][y] == empty_cell) continue;
                int val = input[x][y];
                if (!put(x, y, val)) is_valid = false;
            }
        }
    }
    void init(int siz) {
        set<int> all_numbers;
        for (int v = 1; v <= siz * siz; ++v) all_numbers.insert(v);
        field.assign(siz * siz, vector<int>(siz * siz, UNPUT)),
        nums.assign(siz * siz, vector<vector<int>>(siz * siz, vector<int>(siz * siz, 0))),
        choices.assign(siz * siz, vector<set<int>>(siz * siz, all_numbers)),
        is_valid = true;
    }

    // filed を返す
    const Field& get() {
        return field;
    }
    
    // マス (x, y) に入れられる選択肢を返す
    set<int> find_choices(int x, int y);

    // 数独を出力する
    void print();
    
    // 数独を解く
    void dfs(vector<Field> &res);
    vector<Field> solve();
    
 private:
    // 数独盤面
    int siz;  // default is 3 (9 x 9)
    Field field;

    // nums[x][y][v] ← マス x, y) を含む行・列・ブロックに値 v+1 が何個あるか
    vector<vector<vector<int>>> nums;

    // choices[x][y] ← マス (x, y) に入れることのできる選択肢
    vector<vector<set<int>>> choices;
    
    // 数独の問題として成立しているかどうか
    bool is_valid = true;
    
    // 空きマスのうち、選択肢が最も少ないマスを探す
    bool find_empty(int &x, int &y);

    // マス (x, y) に数値 val を入れることによるマス (x2, y2) への影響
    void put_detail(int x, int y, int val, int x2, int y2);

    // マス (x, y) から数値 val を除去したことによるマス (x2, y2) への影響
    void reset_detail(int x, int y, int val, int x2, int y2);
    
    // マス (x, y) に数値 val を入れる
    bool put(int x, int y, int val);

    // マス (x, y) の数値を削除する
    void reset(int x, int y);

    // 一意に決まるマスを埋めていく
    void process();
};

// マス (x, y) に入れられる選択肢を返す
set<int> Sudoku::find_choices(int x, int y) {
    return choices[x][y];
}

// 空きマスのうち、選択肢が最も少ないマスを探す
bool Sudoku::find_empty(int &x, int &y) {
    const int INF = siz * siz + 1;
    size_t min_num_choices = INF;
    for (int i = 0; i < siz * siz; ++i) {
        for (int j = 0; j < siz * siz; ++j) {
            if (field[i][j] != UNPUT) continue;
            if (min_num_choices > choices[i][j].size()) {
                min_num_choices = choices[i][j].size();
                x = i;
                y = j;
            }
        }
    }
    // 存在しない場合は false
    if (min_num_choices == INF) return false;
    else return true;

}

// マス (x, y) に数値 val を入れることによるマス (x2, y2) への影響
void Sudoku::put_detail(int x, int y, int val, int x2, int y2) {
    if (x == x2 && y == y2) return;
    
    // (x2, y2) にすでに値が入っている場合は何もしない
    if (field[x2][y2] != UNPUT) return;
    
    // それまで (x2, y2) の影響範囲に値 val がなかった場合は選択肢から除く
    if (nums[x2][y2][val - 1] == 0) choices[x2][y2].erase(val);
    
    // nums を更新
    ++nums[x2][y2][val - 1];
}

// マス (x, y) に数値 val を入れる
bool Sudoku::put(int x, int y, int val) {
    // 数値を入れられない場合は false を返す
    if (!find_choices(x, y).count(val)) return false;
    
    // 数値を入れる
    field[x][y] = val;

    // マス (x, y) を含む行・列・ブロックへの影響を更新する
    for (int i = 0; i < siz * siz; ++i) put_detail(x, y, val, x, i);
    for (int i = 0; i < siz * siz; ++i) put_detail(x, y, val, i, y);
    int cx = x / siz * siz + 1, cy = y / siz * siz + 1;
    for (int i = cx - 1; i <= cx + 1; ++i) {
        for (int j = cy - 1; j <= cy + 1; ++j) {
            put_detail(x, y, val, i, j);
        }
    }
    return true;
}

// マス (x, y) から数値 val を除去したことによるマス (x2, y2) への影響
void Sudoku::reset_detail(int x, int y, int val, int x2, int y2) {
    if (x == x2 && y == y2) return;
    
    // (x2, y2) にすでに値が入っている場合は何もしない
    if (field[x2][y2] != UNPUT) return;
    
    // nums を更新
    --nums[x2][y2][val - 1];
    
    // nums が 0 になる場合は選択肢に復活する
    if (nums[x2][y2][val - 1] == 0) choices[x2][y2].insert(val);
}

// マス (x, y) の数値を削除する
void Sudoku::reset(int x, int y) {
    // マス (x, y) を含む行・列・ブロックへの影響を更新する
    int val = field[x][y];
    for (int i = 0; i < siz * siz; ++i) reset_detail(x, y, val, x, i);
    for (int i = 0; i < siz * siz; ++i) reset_detail(x, y, val, i, y);
    int cx = x / siz * siz + 1, cy = y / siz * siz + 1;
    for (int i = cx - 1; i <= cx + 1; ++i) {
        for (int j = cy - 1; j <= cy + 1; ++j) {
            reset_detail(x, y, val, i, j);
        }
    }
    
    // 数値を除去する
    field[x][y] = UNPUT;
}

// 一意に決まるマスを埋めていく
void Sudoku::process() {
    // 数値 1, 2, ..., siz*siz (= default: 9) について順に処理する
    for (int val = 1; val <= siz * siz; ++val) {
        // x 行目について
        for (int x = 0; x < siz * siz; ++x) {
            bool exist = false;
            vector<int> can_enter;
            for (int y = 0; y < siz * siz; ++y) {
                if (field[x][y] == val) exist = true;
                if (field[x][y] == UNPUT
                    && choices[x][y].count(val)) {
                    can_enter.push_back(y);
                }
            }
            // val を入れられるマス目がただ一つならば入れる
            if (!exist && can_enter.size() == 1) {
                int y = can_enter[0];
                put(x, y, val);
            }
        }

        // y 列目について
        for (int y = 0; y < siz * siz; ++y) {
            bool exist = false;
            vector<int> can_enter;
            for (int x = 0; x < siz * siz; ++x) {
                if (field[x][y] == val) exist = true;
                if (field[x][y] == UNPUT
                    && choices[x][y].count(val)) {
                    can_enter.push_back(x);
                }
            }
            // val を入れられるマス目がただ一つならば入れる
            if (!exist && can_enter.size() == 1) {
                int x = can_enter[0];
                put(x, y, val);
            }
        }

        // 各ブロックについて
        for (int bx = 0; bx < siz; ++bx) {
            for (int by = 0; by < siz; ++by) {
                bool exist = false;
                vector<pair<int,int>> can_enter;
                for (int x = bx * siz; x < (bx + 1) * siz; ++x) {
                    for (int y = by * siz; y < (by + 1) * siz; ++y) {
                        if (field[x][y] == val) exist = true;
                        if (field[x][y] == UNPUT
                            && choices[x][y].count(val)) {
                            can_enter.emplace_back(x, y);
                        }
                    }
                }
                // val を入れられるマス目がただ一つならば入れる
                if (!exist && can_enter.size() == 1) {
                    int x = can_enter[0].first;
                    int y = can_enter[0].second;
                    put(x, y, val);
                }
            }
        }
    }
}

// 数独を出力する
void Sudoku::print() {
    for (int x = 0; x < 9; ++x) {
        for (int y = 0; y < 9; ++y) {
            if (field[x][y] == UNPUT) cout << "*";
            else cout << field[x][y];
            cout << " ";
        }
        cout << endl;
    }
}

// 数独を解く
void Sudoku::dfs(vector<Field> &res) {
    // 数独の盤面状態を保持しておく
    const Sudoku board_prev = (*this);
    
    // 一意に自動的に決まるマスを埋める
    process();
    
    // 空きマスの座標を表す
    int x, y;

    // 終端条件を処理し、同時に空きマスを探す
    if (!find_empty(x, y)) {
        // 解に追加
        res.push_back(get());

        // リターンする前に一回元に戻す
        (*this) = board_prev;
        return;
    }
    
    // マス (x, y) に入れられる数値の集合を求める
    const auto &can_use = find_choices(x, y);

    // バックトラッキング
    for (auto val : can_use) {
        put(x, y, val);
        dfs(res);
        reset(x, y);
    }

    // 元に戻す
    (*this) = board_prev;
}

vector<Sudoku::Field> Sudoku::solve() {
    vector<Sudoku::Field> res;
    if (!is_valid) return res;
    dfs(res);
    return res;
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void user_small_test() {
    // 入力
    vector<string> input = {
        "53**7****",
        "6**195***",
        "*98****6*",
        "8***6***3",
        "4**8*3**1",
        "7***2***6",
        "*6****28*",
        "***419**5",
        "****8**79"
    };
    Sudoku sudoku(input);
    
    // 数独を解く
    vector<vector<vector<int>>> res = sudoku.solve();

    // 解を出力する
    if (res.size() == 0) {
        cout << "no solutions." << endl;
    } else {
        if (res.size() > 1) {
            cout << "more than one solutions." << endl;
        }
        cout << "a solution is:" << endl;
        for (int x = 0; x < 9; ++x) {
            for (int y = 0; y < 9; ++y) {
                cout << res[0][x][y] << " ";
            }
            cout << endl;
        }
    }
}

// POJ 2676
void POJ_2676() {
    int T;
    cin >> T;
    
    while (T--) {
        vector<string> input(9);
        for (int i = 0; i < 9; ++i) cin >> input[i];
        Sudoku sudoku(input, '0');
        vector<vector<vector<int>>> res = sudoku.solve();
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) cout << res[0][i][j];
            cout << endl;
        }
    }
}

// ABC 327 C
void ABC_327_C() {
    vector<vector<int>> input(9, vector<int>(9));
    for (int x = 0; x < 9; ++x) for (int y = 0; y < 9; ++y) cin >> input[x][y];
    Sudoku sudoku(input);
    const auto &res = sudoku.solve();
    if (!res.empty()) cout << "Yes" << endl;
    else cout << "No" << endl;
}


int main() {
    user_small_test();
    //POJ_2676();
    //ABC_327_C();
}

