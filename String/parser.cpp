//
// LL1 再帰降下パーサー
//   計算結果を出力するだけでなく、構文解析木を構築する
//
// verified:
//   AOJ 0109 スマート計算機 (for basic expression)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=0109&lang=jp
//
//   AOJ 2613 Unordered Operators (for operator priority)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2613
//
//   AOJ 2401 恒等式 (for unary "not" operator)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2401
//
//   AOJ 1346 Miscalculation (for operator priority)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1346
//
//   TTPC 2023 F - N^a (log N)^b (for unary "log" operator)
//     https://atcoder.jp/contests/ttpc2023/tasks/ttpc2023_f
//
//   AtCoder codeFlyer E - 数式とクエリ (for parse-tree)
//     https://atcoder.jp/contests/bitflyer2018-final-open/tasks/bitflyer2018_final_e
//


/*
 basic settings:
 　・vector<set<string>> OPERATORS
        - 登場する二項演算子と、その優先順位を設定する
 　・function<T(const string&, T, T)> OPERATION
        - 二項演算子に応じた二項演算の結果を出力する処理を設定する
 
 optional settings:
 　・function<T(const string&, int&)> GET_LEAF
        - 構文解析木において葉となる部分のうち、数値以外を処理する分を設定する (第二引数は文字列中の index)
 　・string UNARY_OPERATORS
        - パースする文字列に単項演算子 ("-" や "log" など) が含まれる場合、その関数を表す文字列
 　・function<T(T)> UNARY_OPERATION
        - パースする文字列が関数 OPERATOR_FUNC を含む場合、その関数の中身を設定
 */


#include <bits/stdc++.h>
using namespace std;


// Parser (T must have constructor T(long long))
template<class T> struct Parser {
    /* basic settings */
    vector<vector<string>> OPERATORS;               // binary operator priority
    function<T(const string&, T, T)> OPERATION;     // binary operation

    /* optional settings */
    function<T(const string&, int&)> GET_LEAF;      // leaf parser
    vector<string> UNARY_OPERATORS;                 // unary operator (ex: "-", "log")
    function<T(const string&, T)> UNARY_OPERATION;  // unary operation
    const string OPERATOR_NUMBER = "num";
    
    /* results */
    string S;
    int root;                    // vals[root] is the answer
    vector<T> vals;              // value of each node
    vector<string> ops;          // operator of each node
    vector<int> left, right;     // the index of left-node, right-node
    vector<int> ids;             // the node-index of i-th value
    
    /* constructor */
    Parser() {}
    Parser(const vector<vector<string>> &operators,
           const function<T(const string&, T, T)> operation) {
        init(operators, operation);
    }
    Parser(const vector<vector<string>> &operators,
           const function<T(const string&, T, T)> operation,
           const function<T(const string&, int&)> get_leaf) {
        init(operators, operation, get_leaf);
    }
    Parser(const vector<vector<string>> &operators,
           const function<T(const string&, T, T)> operation,
           const vector<string> &unary_operators,
           const function<T(const string &op, T)> unary_operation) {
        init(operators, operation, unary_operators, unary_operation);
    }
    Parser(const vector<vector<string>> &operators,
           const function<T(const string&, T, T)> operation,
           const function<T(const string&, int&)> get_leaf,
           const vector<string> &unary_operators,
           const function<T(const string &op, T)> unary_operation) {
        init(operators, operation, get_leaf, unary_operators, unary_operation);
    }
    void init(const vector<vector<string>> &operators,
              const function<T(const string&, T, T)> operation) {
        OPERATORS = operators, OPERATION = operation;
    }
    void init(const vector<vector<string>> &operators,
              const function<T(const string&, T, T)> operation,
              const function<T(const string&, int&)> get_leaf) {
        OPERATORS = operators, OPERATION = operation;
        GET_LEAF = get_leaf;
    }
    void init(const vector<vector<string>> &operators,
              const function<T(const string&, T, T)> operation,
              const vector<string> &unary_operators,
              const function<T(const string &op, T)> unary_operation) {
        OPERATORS = operators, OPERATION = operation;
        UNARY_OPERATORS = unary_operators, UNARY_OPERATION = unary_operation;
    }
    void init(const vector<vector<string>> &operators,
              const function<T(const string&, T, T)> operation,
              const function<T(const string&, int&)> get_leaf,
              const vector<string> &unary_operators,
              const function<T(const string &op, T)> unary_operation) {
        OPERATORS = operators, OPERATION = operation;
        GET_LEAF = get_leaf;
        UNARY_OPERATORS = unary_operators, UNARY_OPERATION = unary_operation;
    }
    void clear() {
        S = "";
        vals.clear(), ops.clear(), left.clear(), right.clear(), ids.clear();
    }
    
    /* node generator */
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
    int new_node_with_left(T val, const string &op, int lp) {
        vals.push_back(val);
        ops.push_back(op);
        left.push_back(lp);
        right.push_back(-1);
        ids.push_back(-1);
        return (int)vals.size() - 1;
    }
    // (lp) op (rp)
    int new_node_with_left_right(const string &op, int lp, int rp) {
        vals.push_back(OPERATION(op, vals[lp], vals[rp]));
        ops.push_back(op);
        left.push_back(lp);
        right.push_back(rp);
        ids.push_back(-1);
        return (int)vals.size() - 1;
    }
    
    /* main parser */
    T eval(const string &str, bool default_parse_numeric = true) {
        clear();
        S = str;
        int p = 0, id = 0;
        root = parse(p, id, default_parse_numeric);
        return vals[root];
    }
    int parse(int &p, int &id, bool default_parse_numeric, int depth = 0) {
        assert(p >= 0 && p < (int)S.size());
        string op;
        if (depth < (int)OPERATORS.size()) {
            // recursive descent
            int lp = parse(p, id, default_parse_numeric, depth + 1);
            bool update = false;
            do {
                update = false;
                if (is_in_the_operators(p, op,
                                        OPERATORS[(int)OPERATORS.size() - depth - 1])) {
                    int rp = parse(p, id, default_parse_numeric, depth + 1);
                    lp = new_node_with_left_right(op, lp, rp);
                    update = true;
                }
            } while (p < (int)S.size() && update);
            return lp;
        }
        else if (S[p] == '(') {
            return get_bracket(p, id, default_parse_numeric);
        }
        else if (default_parse_numeric && isdigit(S[p])) {
            return get_number(p, id);
        }
        else if (is_in_the_operators(p, op, UNARY_OPERATORS)) {
            return get_unary_operation(p, id, op, default_parse_numeric);
        }
        else {
            return new_leaf(GET_LEAF(S, p), id);
        }
    }
    bool is_the_operator(int &p, const string &op) {
        if (op != "" && S.substr(p, op.size()) == op) {
            p += op.size();
            return true;
        }
        return false;
    }
    bool is_in_the_operators(int &p, string &res_op, const vector<string> &ops) {
        for (const auto &op : ops) {
            if (is_the_operator(p, op)) {
                res_op = op;
                return true;
            }
        }
        return false;
    }
    int get_number(int &p, int &id) {
        long long val = 0, sign = 1;
        if (S[p] == '-') sign = -1, ++p;
        while (p < (int)S.size() && isdigit(S[p])) {
            val = val * 10 + (int)(S[p++] - '0');
        }
        return new_leaf(T(val * sign), id);
    }
    int get_bracket(int &p, int &id, bool default_parse_numeric) {
        ++p;  // skip "("
        int lp = parse(p, id, default_parse_numeric, 0);
        ++p;  // skip ")"
        return lp;
    }
    int get_unary_operation(int &p, int &id, const string &op,
                            bool default_parse_numeric) {
        int lp;
        if (S[p] == '(') lp = get_bracket(p, id, default_parse_numeric);
        else lp = parse(p, id, default_parse_numeric, (int)OPERATORS.size());
        return new_node_with_left(UNARY_OPERATION(op, vals[lp]), op, lp);
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

/* user small test */
void small_test() {
    // set parser
    vector<vector<string>> operators = {{"*", "/"}, {"+", "-"}};
    auto operation = [&](const string &op, long long a, long long b) -> long long {
        if (op == "+") return a + b;
        else if (op == "-") return a - b;
        else if (op == "*") return a * b;
        else if (op == "/") return a / b;
        else return 0;
    };
    vector<string> unary_operators = {"-"};
    auto unary_operation = [&](const string &op, long long a) -> long long {
        return -a;
    };
    Parser<long long> ps(operators, operation, unary_operators, unary_operation);
    
    vector<string> tests = {
        "-(-3)",
        "-(-9+3*(-4))",
        "-5-(-3)",
        "-4-99"
    };
    for (const auto &S : tests) {
        cout << S << " = " << ps.eval(S) << endl;
    }
}


/* AOJ 0109 */
void AOJ_0109() {
    vector<vector<string>> operators = {{"*", "/"}, {"+", "-"}};
    auto operation = [&](const string &op, long long a, long long b) -> long long {
        if (op == "+") return a + b;
        else if (op == "-") return a - b;
        else if (op == "*") return a * b;
        else if (op == "/") return a / b;
        else return 0;
    };
    Parser<long long> ps(operators, operation);
    
    int N;
    cin >> N;
    for (int i = 0; i < N; ++i) {
        string S;
        cin >> S;
        S.pop_back();
        cout << ps.eval(S) << endl;
    }
}


/* AOJ 2613 */
void AOJ_2613() {
    string S;
    cin >> S;
    
    // 二項演算を設定
    auto operation = [&](const string &op, long long l, long long r) -> long long {
        if (op == "+") return l + r;
        else if (op == "-") return l - r;
        else return l * r;
    };

    // 演算子の優先順位をすべて試す
    long long res = -LONG_MAX;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                vector<vector<string>> operators(3);
                operators[i].push_back("+");
                operators[j].push_back("-");
                operators[k].push_back("*");
                Parser<long long> ps(operators, operation);
                res = max(res, ps.eval(S));
            }
        }
    }
    cout << res << endl;
}
    

/* AOJ 2401 */
void AOJ_2401() {
    // parser settings //
    vector<vector<string>> ops = {{"+", "*", "->"}};
    auto operation = [&](const string &op, bool l, bool r) -> bool {
        if (op == "+") return (l || r);
        else if (op == "*") return (l && r);
        else if (op == "->") return (!l || r);
        else return true;
    };
    auto get_leaf = [&](const string &S, int &p) -> bool {
        return (S[p++] == 'T' ? true : false);
    };
    vector<string> unary_ops = {"-"};
    auto unary_operation = [&](const string &op, bool l) -> bool {
        if (op == "-") return !l;
        else return l;
    };
    Parser<bool> ps1(ops, operation, get_leaf, unary_ops, unary_operation);
    Parser<bool> ps2(ops, operation, get_leaf, unary_ops, unary_operation);
    
    string S;
    while (cin >> S) {
        if (S == "#") break;
        for (int i = 0; i < S.size(); ++i) if (S[i] == '=') S[i] = ' ';
        
        bool res = true;
        for (int bit = 0; bit < (1<<11); ++bit) {
            string SS = S;
            for (int i = 0; i < SS.size(); ++i) {
                if (SS[i] >= 'a' && SS[i] <= 'k') {
                    int id = SS[i] - 'a';
                    if (bit & (1<<id)) SS[i] = 'T';
                    else SS[i] = 'F';
                }
            }
            istringstream si(SS);
            string s1, s2;
            si >> s1 >> s2;
            bool r1 = ps1.eval(s1, false), r2 = ps2.eval(s2, false);
            if (r1 != r2) {
                res = false;
            }
        }
        if (res) cout << "YES" << endl;
        else cout << "NO" << endl;
    }
}


/* AOJ 1346 */
void AOJ_1346() {
    // multiplication-rool parser, left-to-right parser
    vector<vector<string>> operator_mul = {{"*"}, {"+"}};
    vector<vector<string>> operator_left_right = {{"+", "*"}};
    auto operation = [&](string op, long long a, long long b) -> long long {
        if (op == "+") return a + b;
        else return a * b;
    };
    Parser<long long> ps_mul(operator_mul, operation);
    Parser<long long> ps_left_right(operator_left_right, operation);
    
    string S;
    long long ans;
    cin >> S >> ans;
    
    bool mul = false, left_right = false;
    if (ps_mul.eval(S) == ans) mul = true;
    if (ps_left_right.eval(S) == ans) left_right = true;
    
    if (mul && left_right) cout << "U" << endl;
    else if (mul) cout << "M" << endl;
    else if (left_right) cout << "L" << endl;
    else cout << "I" << endl;
}


/* TTPC 2023 F */
struct Node {
    long long n, l, val;
    bool eps;
    Node(long long n = 0, long long l = 0, long long e = 0, long long v = 0)
        : n(n), l(l), eps(e), val(v) {}
    friend ostream& operator << (ostream &s, Node node) {
        return s << "(" << node.n << ", " << node.l << ", " << node.eps
                << ", " << node.val << ")";
    }
};

void TTPC_2023_F() {
    // operators
    vector<vector<string>> operators = {{"^"}, {"*"}, {"+"}};
    
    // define binary operation
    auto operation = [&](const string &op, const Node &l, const Node &r) -> Node {
        if (op == "+") {
            vector<long long> vl({l.n, l.l, (long long)l.eps, l.val});
            vector<long long> vr({r.n, r.l, (long long)r.eps, r.val});
            return (vl > vr ? l : r);
        }
        else if (op == "*") return Node(l.n + r.n, l.l + r.l, l.eps | r.eps, 1);
        else if (op == "^") return Node(l.n * r.val, l.l * r.val, l.eps, 1);
        else return Node();
    };
    
    // define number process
    auto get_leaf = [&](const string &S, int &p) -> Node {
        if (S[p] == 'N') {
            ++p;
            return Node(1, 0, false, 1);
        } else {
            long long val = 0;
            do {
                val = val * 10 + (int)(S[p++] - '0');
            } while (S[p] >= '0' && S[p] <= '9');
            return Node(0, 0, false, val);
        }
    };
    
    // define function (log)
    vector<string> unary_operators = {"log"};
    auto unary_operation = [&](const string &op, const Node &node) -> Node {
        if (op == "log") {
            if (node.n >= 1) return Node(0, 1, false, 1);
            else return Node(0, 0, true, 0);
        } else return Node();
    };
    
    // 入力
    string S;
    cin >> S;
    Parser<Node> ps(operators, operation, get_leaf, unary_operators, unary_operation);
    auto res = ps.eval(S, false);
    //ps.dump();

    if (res.eps) ++res.l;
    cout << res.n << " " << res.l << endl;
}


/* codeFloyer E */
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() : val(0) { }
    constexpr Fp(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr long long get() const { return val; }
    constexpr int get_mod() const { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp(0) - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        val = val * u % MOD;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp pow(long long n) const {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        Fp res(1), div(*this);
        return res / div;
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) const {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) const {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) {
        return os << x.val;
    }
    friend constexpr Fp<MOD> pow(const Fp<MOD> &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> inv(const Fp<MOD> &r) {
        return r.inv();
    }
};
 
void codeFloyer_E() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    
    // 入力
    string S;
    int Q;
    cin >> S >> Q;
    int N = 0;
    for (int i = 0; i < S.size(); ++i) if (S[i] == 'a') ++N;
    vector<long long> given(N);
    for (int i = 0; i < N; ++i) cin >> given[i];
    
    // set Parser
    vector<vector<string>> operators = {{"*"}, {"+", "-"}};
    auto operation = [&](const string &op, mint a, mint b) -> mint {
        if (op == "+") return a + b;
        else if (op == "-") return a - b;
        else if (op == "*") return a * b;
        else return 0;
    };
    int id = 0;
    auto get_leaf = [&](const string &S, int &p) -> mint {
        if (S[p] == 'a') {
            ++p;
            return given[id++];
        } else return 0;
    };
    Parser<mint> ps(operators, operation, get_leaf);
    mint base = ps.eval(S);
    
    // DP
    vector<mint> dp(N);
    auto rec = [&](auto self, int v, mint w) -> void {
        if (ps.ops[v] == ps.OPERATOR_NUMBER) dp[ps.ids[v]] = w;
        else if (ps.ops[v] == "+") {
            self(self, ps.left[v], w);
            self(self, ps.right[v], w);
        }
        else if (ps.ops[v] == "-") {
            self(self, ps.left[v], w);
            self(self, ps.right[v], -w);
        }
        else if (ps.ops[v] == "*") {
            self(self, ps.left[v], w * ps.vals[ps.right[v]]);
            self(self, ps.right[v], w * ps.vals[ps.left[v]]);
        }
    };
    rec(rec, ps.root, mint(1));
    
    // 出力
    for (int q = 0; q < Q; ++q) {
        long long id, x;
        cin >> id >> x;
        --id;
        mint res = base + dp[id] * (x - given[id]);
        cout << res << endl;
    }
}


int main() {
    small_test();
    //AOJ_0109();
    //AOJ_2613();
    //AOJ_2401();
    //AOJ_1346();
    //TTPC_2023_F();
    //codeFloyer_E();
}



