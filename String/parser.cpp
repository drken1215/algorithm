//
// LL1 再帰降下パーサー
//
// verified:
//   TTPC 2023 F - N^a (log N)^b
//     https://atcoder.jp/contests/ttpc2023/tasks/ttpc2023_f
//


/*
 setting:
 　・vector<set<char>> OPERATORS - 登場する演算子と、その優先順位を設定
 　・function<T(char, T, T)> OP - 演算子に応じた二項演算の結果を出力する処理を設定
 　・function<T(const string&, int&)> GETNUM - パースする文字列中の数字部分を抽出する処理を設定
                                              (第二引数は文字列中の index を表す)
 　・function<T(T)> FUNC - パースする文字列が関数 "f" を含む場合、その関数の中身を設定
 */


#include <bits/stdc++.h>
using namespace std;


// Parser
template<class T> struct Parser {
    /* settings */
    vector<set<char>> OPERATORS;              // operator priority
    function<T(char, T, T)> OP;               // binary operation
    function<T(const string&, int&)> GETNUM;  // number parser
    function<T(T)> FUNC;                      // inner function in the string
    
    const char OPERATOR_NUMBER = 'n';
    const char OPERATOR_FUNC = 'f';
    
    /* results */
    int root;                    // vals[root] is the answer
    vector<T> vals;              // value of each node
    vector<char> ops;            // operator of each node ('a' means leaf values)
    vector<int> left, right;     // the index of left-node, right-node
    vector<int> ids;             // the node-index of i-th value
    
    /* constructor */
    Parser() {}
    Parser(const vector<set<char>> &operators,
           const function<T(char, T, T)> op,
           const function<T(const string&, int&)> getnum) {
        init(operators, op, getnum);
    }
    Parser(const vector<set<char>> &operators,
           function<T(char, T, T)> op,
           const function<T(const string&, int&)> getnum,
           const function<T(T)> func) {
        init(operators, op, getnum, func);
    }
    void init(const vector<set<char>> &operators,
              const function<T(char, T, T)> op,
              const function<T(const string&, int&)> getnum) {
        OPERATORS = operators;
        OP = op;
        GETNUM = getnum;
    }
    void init(const vector<set<char>> &operators,
              const function<T(char, T, T)> op,
              const function<T(const string&, int&)> getnum,
              const function<T(T)> func) {
        OPERATORS = operators;
        OP = op;
        GETNUM = getnum;
        FUNC = func;
    }
    void clear() {
        vals.clear(), ops.clear(), left.clear(), right.clear(), ids.clear();
    }
    T parse(const string &S) {
        clear();
        int p = 0, id = 0;
        root = parse(S, p, id);
        return vals[root];
    }
    
    /* main parser */
    int parse(const string &S, int &p, int &id, int depth = 0) {
        if (depth < (int)OPERATORS.size()) {
            int lp = parse(S, p, id, depth + 1);
            while (p < (int)S.size() &&
                   OPERATORS[(int)OPERATORS.size() - depth - 1].count(S[p])) {
                char op = S[p++];
                int rp = parse(S, p, id, depth + 1);
                lp = new_node_with_left_right(op, lp, rp);
            }
            return lp;
        } else if (S[p] == '(') {
            ++p;                        // skip "("
            int lp = parse(S, p, id, 0);
            ++p;                        // skip ")"
            return lp;
        } else if (S[p] == OPERATOR_FUNC) {
            p += 2;                     // skip "f("
            int lp = parse(S, p, id, 0);
            ++p;                        // skip ")"
            return new_node_with_left(lp);
        } else {
            return new_leaf(GETNUM(S, p), id);
        }
    }
    
    /* generate node */
    // value
    int new_leaf(T val, int &id) {
        vals.push_back(val);
        ops.push_back(OPERATOR_NUMBER);
        left.push_back(-1);
        right.push_back(-1);
        ids.push_back(id++);
        return (int)vals.size() - 1;
    }
    // FUNC(lp)
    int new_node_with_left(int lp) {
        vals.push_back(FUNC(vals[lp]));
        ops.push_back(OPERATOR_FUNC);
        left.push_back(lp);
        right.push_back(-1);
        ids.push_back(-1);
        return (int)vals.size() - 1;
    }
    // (lp) op (rp)
    int new_node_with_left_right(char op, int lp, int rp) {
        vals.push_back(OP(op, vals[lp], vals[rp]));
        ops.push_back(op);
        left.push_back(lp);
        right.push_back(rp);
        ids.push_back(-1);
        return (int)vals.size() - 1;
    }
    
    /* dump */
    void dump() {
        dump_rec(root);
    }
    void dump_rec(int v, string space = "") {
        cout << space << v << ": (" << ops[v] << ", " << vals[v]
             << ") -> left: " << left[v] << ", right: " << right[v] << endl;
        if (left[v] != -1) dump_rec(left[v], space + "  ");
        if (right[v] != -1) dump_rec(right[v], space + "  ");
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

/* TTPC 2023 F */
struct Node {
    long long n, l, val;
    bool eps;
    Node() {}
    Node(long long n, long long l, long long e, long long v)
        : n(n), l(l), eps(e), val(v) {}
};

ostream& operator << (ostream &s, Node node) {
    return s << "(" << node.n << ", " << node.l << ", " << node.eps << ", " << node.val << ")";
}

void TTPC_2023_F() {
    // operators
    vector<set<char>> ops = {
        set<char>({'^'}),
        set<char>({'*'}),
        set<char>({'+'})
    };
    
    // define binary operation
    auto op = [&](char op, const Node &l, const Node &r) -> Node {
        if (op == '+') {
            vector<long long> vl({l.n, l.l, (long long)l.eps, l.val});
            vector<long long> vr({r.n, r.l, (long long)r.eps, r.val});
            return (vl > vr ? l : r);
        }
        else if (op == '*') return Node(l.n + r.n, l.l + r.l, l.eps | r.eps, 1);
        else if (op == '^') return Node(l.n * r.val, l.l * r.val, l.eps, 1);
        else return Node();
    };
    
    // define number process
    auto getnum = [&](const string &S, int &p) -> Node {
        Node res;
        if (S[p] == 'N') {
            ++p;
            return Node(1, 0, false, 1);
        } else {
            long long val = 0;
            while (S[p] >= '0' && S[p] <= '9') {
                val *= 10;
                val += (int)(S[p] - '0');
                ++p;
            }
            return Node(0, 0, false, val);
        }
    };
    
    // define function (log)
    auto func = [&](const Node &node) -> Node {
        if (node.n >= 1) return Node(0, 1, false, 1);
        else return Node(0, 0, true, 0);
    };
    
    // 入力, "log()" -> "f()"
    string S;
    cin >> S;
    string T = "";
    for (int i = 0; i < S.size(); ++i) {
        if (S[i] == 'l') T += "f", i += 2;
        else T += S[i];
    }
    
    Parser<Node> ps(ops, op, getnum, func);
    auto res = ps.parse(T);
    //ps.dump();

    if (res.eps) ++res.l;
    cout << res.n << " " << res.l << endl;
}


int main() {
    TTPC_2023_F();
}

