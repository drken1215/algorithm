//
// 多倍長整数
//
// verified
//   ARC 007 D - 破れた宿題
//     https://atcoder.jp/contests/arc007/tasks/arc007_4
//


#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;


const int DEFAULT_SIZE = 125;
struct bint : vector<long long> {
    static const long long BASE = 100000000;
    static const int BASE_DIGIT = 8;
    int sign;

    // constructor
    bint(long long num = 0) : vector<long long>(DEFAULT_SIZE, 0), sign(1) {
        if (num < 0) sign = -1, num = -num;
        (*this)[0] = num; 
        this->normalize();
    }
    bint(int size, long long num) : vector<long long>(size, num), sign(1) {}
    bint& normalize() {
        long long c = 0;
        bool exist = false;
        for (int i = 0;; ++i) {
            if (i >= this->size()) this->push_back(0);
            if ((*this)[i] < 0 && i+1 >= this->size()) this->push_back(0);
            while ((*this)[i] < 0) {
                (*this)[i+1] -= 1;
                (*this)[i] += BASE;
            }
            long long a = (*this)[i] + c;
            (*this)[i] = a % BASE;
            if ((*this)[i]) exist = true;
            c = a / BASE;
            if (c == 0 && i == this->size()-1) break;
        }
        if (!exist) sign = 1;
        return (*this);
    }
    friend bint abs(const bint &x) {
        bint z = x;
        if (z.sign == -1) z.sign = 1;
        return z;
    }

    // operation
    bint operator - () const {
        bint res = *this;
        bool allzero = true;
        for (int i = 0; i < this->size(); ++i) {
            if (res[i] != 0) {
                allzero = false;
                break;
            }
        }
        if (!allzero) res.sign = -res.sign; 
        return res; 
    }
    bint& operator += (const bint& r) {
        while (size() < r.size()) this->emplace_back(0);
        if (sign == r.sign) {
            for (int i = 0; i < r.size(); ++i) (*this)[i] += r[i];
        }
        else {
            if (sign == 1 && abs(*this) < abs(r)) sign = -1;
            else if (sign == -1 && abs(*this) <= abs(r)) sign = 1;
            if (abs(*this) >= abs(r)) {
                for (int i = 0; i < r.size(); ++i) (*this)[i] -= r[i];
            }
            else {
                for (int i = 0; i < size(); ++i) (*this)[i] = -(*this)[i];
                for (int i = 0; i < r.size(); ++i) (*this)[i] += r[i];
            }
        }
        return this->normalize();
    }
    bint& operator -= (const bint& r) {
        while (size() < r.size()) this->emplace_back(0);
        if (sign == -r.sign) {
            for (int i = 0; i < r.size(); ++i) (*this)[i] += r[i];
        }
        else {
            if (sign == 1 && abs(*this) < abs(r)) sign = -1;
            else if (sign == -1 && abs(*this) <= abs(r)) sign = 1;
            if (abs(*this) >= abs(r)) {
                for (int i = 0; i < r.size(); ++i) (*this)[i] -= r[i];
            }
            else {
                for (int i = 0; i < size(); ++i) (*this)[i] = -(*this)[i];
                for (int i = 0; i < r.size(); ++i) (*this)[i] += r[i];
            }
        }
        return this->normalize();
    }
    bint& operator *= (long long r) {
        if ( (sign == 1 && r >= 0) || (sign == -1 && r < 0) ) sign = 1;
        else sign = -1;
        if (r < 0) r = -r;
        for (int i = 0; i < size(); ++i) (*this)[i] *= r;
        return this->normalize();
    }
    bint& operator *= (const bint& r) {
        int tx = (int)size()-1, ty = (int)r.size()-1; 
        for (tx = size()-1; tx >= 0; --tx) if ((*this)[tx] > 0) break;
        for (ty = r.size()-1; ty >= 0; --ty) if (r[ty] > 0) break;
        bint res(0);
        res.resize(tx+ty+2);
        if (sign == r.sign) res.sign = 1;
        else res.sign = -1;
        for (int i = 0; i <= tx; ++i) {
            for (int j = 0; j <= ty && i+j < (int)res.size()-1; ++j) {
                long long val = (*this)[i] * r[j] + res[i+j];
                res[i+j+1] += val / bint::BASE;
                res[i+j] = val % bint::BASE;
            }
        }
        return (*this) = res.normalize();
    }
    friend bint pow(const bint& a, long long n) {
        bint res(1), b = a;
        while (n > 0) {
            if (n & 1) res = res * b;
            b = b * b;
            n >>= 1;
        }
        return res;
    }
    bint operator + (const bint& r) const { return bint(*this) += r; }
    bint operator - (const bint& r) const { return bint(*this) -= r; }
    bint operator * (long long r) const { return bint(*this) *= r; }
    bint operator * (const bint& r) const { return bint(*this) *= r; }

    // divide
    bint& operator /= (long long r) {
        if (r < 0) sign *= -1, r = -r;
        long long c = 0, t = 0;
        for (int i = (int)size()-1; i >= 0; --i) {
            t = bint::BASE * c + (*this)[i];
            (*this)[i] = t / r;
            c = t % r;
        }
        this->normalize();
        return (*this);
    }
    long long operator %= (long long r) {
        if (r < 0) sign *= -1, r = -r;
        long long c = 0, t = 0;
        for (int i = (int)size()-1; i >= 0; --i) {
            t = bint::BASE * c + (*this)[i];
            (*this)[i] = t / r;
            c = t % r;
        }
        return c;
    }
    bint operator / (long long r) const {
        return bint(*this) /= r;
    }
    long long operator % (long long r) const {
        return bint(*this) %= r;
    }
    friend pair<bint, bint> divmod(const bint &a, const bint &r) {
        bint zero = 0, s = 0, t = 0;
        if (abs(a) < abs(r)) return {zero, a};
        bint ar = abs(r);
        s.resize((int)a.size()), t.resize((int)r.size());
        int tx = (int)a.size()-1;
        for (;tx >= 0; --tx) if (a[tx] > 0) break;
        for (int i = tx; i >= 0; --i) {
            t = t * bint::BASE + a[i];
            long long lo = 0, hi = bint::BASE;
            if (t >= ar) {
                while (hi - lo > 1) {
                    int mid = (hi + lo) / 2;
                    if (ar * mid > t) hi = mid;
                    else lo = mid;
                }
                t -= ar * lo;
            }
            s[i] = lo;
        }
        if (a.sign == r.sign) s.sign = 1, t.sign = 1;
        else s.sign = -1, t.sign = 1;
        return make_pair(s.normalize(), t.normalize());
    }
    bint operator / (const bint& r) const {
        return divmod((*this), r).first;
    }
    bint operator % (const bint& r) const {
        return divmod((*this), r).second;
    }
    bint& operator /= (const bint& r) { return (*this) = (*this) / r; }
    bint& operator %= (const bint& r) { return (*this) = (*this) % r; }

    // equality
    friend bool operator < (const bint &x, const bint& y) {
        if (x.sign < y.sign) return true;
        else if (x.sign > y.sign) return false;
        else {
            int tx = (int)x.size()-1, ty = (int)y.size()-1; 
            for (tx = x.size()-1; tx >= 0; --tx) if (x[tx] > 0) break;
            for (ty = y.size()-1; ty >= 0; --ty) if (y[ty] > 0) break;
            if (tx < ty) return true;
            else if (tx > ty) return false;
            else if (x.sign == 1) {
                for (int i = tx; i >= 0; --i)
                    if (x[i] != y[i]) return x[i] < y[i];
                return false;
            }
            else {
                for (int i = tx; i >= 0; --i)
                    if (x[i] != y[i]) return x[i] > y[i];
                return false;
            }
        }
    }
    friend bool operator > (const bint& x, const bint& y) { return y < x; }
    friend bool operator <= (const bint& x, const bint& y) { return !(x > y); }
    friend bool operator >= (const bint& x, const bint& y) { return !(x < y); }
    friend bool operator == (const bint &x, const bint& y) {
        if (x.sign != y.sign) return false;
        int tx = (int)x.size()-1, ty = (int)y.size()-1; 
        for (tx = x.size()-1; tx >= 0; --tx) if (x[tx] > 0) break;
        for (ty = y.size()-1; ty >= 0; --ty) if (y[ty] > 0) break;
        if (tx != ty) return false;
        for (int i = tx; i >= 0; --i)
            if (x[i] != y[i]) return false;
        return true;
    }
    friend bool operator != (const bint& x, const bint& y) { return !(x == y); }
};

bint toBint(const string &is) {
    string s = is;
    if (s[0] == '-') s = s.substr(1);
    while (s.size() % 8 != 0) s = "0" + s;
    int N = (int)s.size();
    bint res(N/8, 0);
    for (int i = 0; i < (int)s.size(); ++i) {
        res[(N-i-1)/8] *= 10;
        res[(N-i-1)/8] += s[i] - '0';
    }
    if (is[0] == '-') res.sign = -1;
    return res;
}

string toStr(const bint &r) {
    stringstream ss;
    if (r.sign == -1) ss << '-';
    int d = (int)r.size()-1; 
    for (; d >= 0; --d) if (r[d] > 0) break;
    if (d == -1) ss << 0;
    else ss << r[d];
    for (int i = d-1; i >= 0; --i) {
        ss.width(bint::BASE_DIGIT);
        ss.fill('0');
        ss << r[i];
    }
    return ss.str();
}

istream &operator >> (istream &is, bint &x) {
    string s; is >> s;
    x = toBint(s);
    return is;
}

ostream &operator << (ostream &os, const bint &x) {
    if (x.sign == -1) os << '-';
    int d = x.size()-1; 
    for (d = x.size()-1; d >= 0; --d) if (x[d] > 0) break;
    if (d == -1) os << 0;
    else os << x[d];
    for (int i = d-1; i >= 0; --i) {
        os.width(bint::BASE_DIGIT);
        os.fill('0');
        os << x[i];
    }
    return os;
}



// 初項が syoko, 第二項を niko とすることが可能かどうか
// S の先頭は niko の先頭から始まる
bool isValid(const string &S, const bint &syoko, const bint &niko) {
    if (niko <= syoko) return false;
    int N = (int)S.size();
    string sniko = toStr(niko);
    int si = 0;
    for (int i = 0; i < sniko.size() && si < N; ++i) {
        if (sniko[i] != S[si++]) return false;
    }
    if (si == N) return true;

    bint prev = syoko, cur = niko;
    while (si < N) {
        if (S[si] == '0') return false;
        bint next = cur * 2 - prev;
        string snext = toStr(next);
        for (int i = 0; i < snext.size() && si < N; ++i) {
            if (snext[i] != S[si++]) return false;
        }
        if (si < N) prev = cur, cur = next;
    }
    return true;
}

// 試す
void judge(const string &S, const bint& syoko, const bint& niko_koho, bint& niko) {
    if (!isValid(S, syoko, niko_koho)) return;
    if (niko == 0) niko = niko_koho;
    else if (niko_koho < niko) niko = niko_koho;
}

// 解く
void solve(string S) {
    // 初項
    if (S[0] == '0') S = "1" + S;
    int zeronum = 0;
    for (int i = 1; i < S.size(); ++i) {
        if (S[i] != '0') break;
        ++zeronum;
    }
    bint syoko = toBint(S.substr(0, zeronum + 1));
    S = S.substr(zeronum+1);
    int N = (int)S.size();
    if (N == 0) {
        cout << syoko << " " << 1 << endl;
        return;
    }

    bint niko = 0;
    for (int i = 1; i <= N; ++i) {
        bint niko_koho = toBint(S.substr(0, i));
        if (niko_koho <= syoko) continue;
        judge(S, syoko, niko_koho, niko);
    }

    // +1 が valid か
    judge(S, syoko, syoko + 1, niko);

    // 0 を繋げていく
    bint niko_koho = toBint(S);
    while (niko_koho <= syoko) niko_koho *= 10;
    judge(S, syoko, niko_koho, niko);

    // 出力
    cout << syoko << " " << niko - syoko << endl;
}

int main() {
    string S;
    while (cin >> S) solve(S);
}
