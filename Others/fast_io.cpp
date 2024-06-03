//
// Fast IO
//
// references:
//   maspy: [Library Checker] Many A + B
//     https://maspypy.com/library-checker-many-a-b
//
// verified:
//   Yosupo Library Checker - Many A + B
//     https://judge.yosupo.jp/problem/many_aplusb
//


#include <bits/stdc++.h>
using namespace std;


namespace mylib {
struct FastInput {
private:
    static const unsigned int INQSZ = 1 << 17;
    int p = 0;
    int q = 0;
    static char Q[INQSZ + 1];
    void seek1KiB() {
        if (q < p + 1024) {
            if (p <= q - p) memcpy(Q, Q + p, q - p);
            else memmove(Q, Q + p, q - p);
            q -= p;
            p = 0;
            q += fread(Q + q, 1, INQSZ - q, stdin);
            if (q <= 1024) memset(Q + q, 0, 1025 - q);
        }
    }
    
public:
    char seekChar() {
        if (p == q) {
            q = fread(Q, 1, INQSZ, stdin);
            Q[q] = '\0';
            p = 0;
        }
        return Q[p];
    }
    void skipSpace() {
        bool f = true;
        while (f) {
            seek1KiB();
            for (int i = 0; i < 1024; ++i) if (!isspace(Q[p++])) {
                f = false;
                break;
            }
        }
        --p;
    }
    unsigned int nextU32() {
        skipSpace();
        seek1KiB();
        unsigned int buf = 0;
        while (true) {
            char tmp = Q[p];
            if ('9' < tmp || tmp < '0') break;
            buf = buf * 10 + (tmp - '0');
            ++p;
        }
        return buf;
    }
    int nextI32() {
        skipSpace();
        if (Q[p] == '-') {
            ++p;
            return (int)(1 + ~nextU32());
        }
        return (int)nextU32();
    }
    unsigned long long nextU64() {
        skipSpace();
        seek1KiB();
        unsigned long long buf = 0;
        while (true) {
            char tmp = Q[p];
            if ('9' < tmp || tmp < '0') break;
            buf = buf * 10 + (tmp - '0');
            ++p;
        }
        return buf;
    }
    long long nextI64() {
        skipSpace();
        seek1KiB();
        if (seekChar() == '-') {
            ++p;
            return (long long)(-nextU64());
        }
        return (long long)nextU64();
    }
    char nextChar() {
        skipSpace();
        char buf = seekChar();
        ++p;
        return buf;
    }
    string nextToken() {
        skipSpace();
        string buf;
        while (true) {
            char ch = seekChar();
            if (isspace(ch) || ch == '\0') break;
            buf.push_back(ch);
            ++p;
        }
        return buf;
    }
    FastInput& operator >> (unsigned int &dest) { dest = nextU32(); return *this; }
    FastInput& operator >> (int &dest) { dest = nextI32(); return *this; }
    FastInput& operator >> (unsigned long &dest) { dest = nextU64(); return *this; }
    FastInput& operator >> (long &dest) { dest = nextI64(); return *this; }
    FastInput& operator >> (unsigned long long &dest) { dest = nextU64(); return *this; }
    FastInput& operator >> (long long &dest) { dest = nextI64(); return *this; }
    FastInput& operator >> (string &dest) { dest = nextToken(); return *this; }
    FastInput& operator >> (char &dest) { dest = nextChar(); return *this; }
} cin;

struct FastOutputTable {
    char LZ[1000][4] = {};
    char NLZ[1000][4] = {};
    constexpr FastOutputTable() {
        for (uint_fast32_t d = 0; d < 1000; ++d) {
            LZ[d][0] = ('0' + d / 100 % 10);
            LZ[d][1] = ('0' + d /  10 % 10);
            LZ[d][2] = ('0' + d /   1 % 10);
            LZ[d][3] = '\0';
        }
        for (uint_fast32_t d = 0; d < 1000; ++d) {
            uint_fast32_t i = 0;
            if (d >= 100) NLZ[d][i++] = ('0' + d / 100 % 10);
            if (d >=  10) NLZ[d][i++] = ('0' + d /  10 % 10);
            if (d >=   1) NLZ[d][i++] = ('0' + d /   1 % 10);
            NLZ[d][i++] = '\0';
        }
    }
};

struct FastOutput {
private:
    static const int OUTQSZ = 1 << 17;
    static char Q[OUTQSZ];
    static constexpr FastOutputTable TB = FastOutputTable();
    int p = 0;
    static constexpr uint_fast32_t P10(uint_fast32_t d) {
        return d ? P10(d - 1) * 10 : 1;
    }
    static constexpr uint_fast64_t P10L(uint_fast32_t d) {
        return d ? P10L(d - 1) * 10 : 1;
    }
    template<class T, class U> static void Fil(T& m, U& l, U x) {
        m = l/x;
        l -= m*x;
    }
    void next_dig9(uint_fast32_t x) {
        uint_fast32_t y;
        Fil(y, x, P10(6));
        nextCstr(TB.LZ[y]);
        Fil(y, x, P10(3));
        nextCstr(TB.LZ[y]); nextCstr(TB.LZ[x]);
    }
    
public:
    void nextChar(char c) {
        Q[p++] = c;
        if (p == OUTQSZ) {
            fwrite(Q, p, 1, stdout);
            p = 0;
        }
    }
    void nextEoln() {
        nextChar('\n');
    }
    void nextCstr(const char *s) {
        while (*s) nextChar(*(s++));
    }
    void nextU32(uint_fast32_t x) {
        uint_fast32_t y = 0;
        if (x >= P10(9)) {
            Fil(y, x, P10(9));
            nextCstr(TB.NLZ[y]); next_dig9(x);
        } else if (x >= P10(6)) {
            Fil(y, x, P10(6));
            nextCstr(TB.NLZ[y]);
            Fil(y, x, P10(3));
            nextCstr(TB.LZ[y]); nextCstr(TB.LZ[x]);
        } else if (x >= P10(3)) {
            Fil(y, x, P10(3));
            nextCstr(TB.NLZ[y]); nextCstr(TB.LZ[x]);
        }
        else if (x >= 1) nextCstr(TB.NLZ[x]);
        else nextChar('0');
    }
    void nextI32(int_fast32_t x) {
        if (x >= 0) nextU32(x);
        else {
            nextChar('-');
            nextU32((uint_fast32_t)-x);
        }
    }
    void nextU64(uint_fast64_t x) {
        uint_fast32_t y = 0;
        if (x >= P10L(18)) {
            Fil(y, x, P10L(18));
            nextU32(y);
            Fil(y, x, P10L(9));
            next_dig9(y);
            next_dig9(x);
        } else if (x >= P10L(9)) {
            Fil(y, x, P10L(9));
            nextU32(y);
            next_dig9(x);
        } else nextU32(x);
    }
    void nextI64(int_fast64_t x) {
        if (x >= 0) nextU64(x);
        else {
            nextChar('-');
            nextU64((uint_fast64_t)-x);
        }
    }
    void writeToFile(bool flush = false) {
        fwrite(Q, p, 1, stdout);
        if (flush) fflush(stdout);
        p = 0;
    }
    FastOutput() { Q[0] = 0; }
    ~FastOutput() { writeToFile(); }
    FastOutput& operator << (unsigned int tg) { nextU32(tg); return *this; }
    FastOutput& operator << (unsigned long tg) { nextU64(tg); return *this; }
    FastOutput& operator << (unsigned long long tg) { nextU64(tg); return *this; }
    FastOutput& operator << (int tg) { nextI32(tg); return *this; }
    FastOutput& operator << (long tg) { nextI64(tg); return *this; }
    FastOutput& operator << (long long tg) { nextI64(tg); return *this; }
    FastOutput& operator << (const string& tg) { nextCstr(tg.c_str()); return *this; }
    FastOutput& operator << (const char* tg) { nextCstr(tg); return *this; }
    FastOutput& operator << (char tg) { nextChar(tg); return *this; }
} cout;

char FastInput::Q[INQSZ + 1];
char FastOutput::Q[OUTQSZ];
} // namespace mylib



//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Many A + B
void Yosupo_A_puls_B() {
    using mylib::cin, mylib::cout;
    
    int T;
    cin >> T;
    for (int t = 0; t < T; ++t) {
        unsigned long long a, b;
        cin >> a >> b;
        unsigned long long ans = a + b;
        cout << ans << '\n';
    }
}


int main() {
    Yosupo_A_puls_B();
}

