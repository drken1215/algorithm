//
// 任意 mod の畳み込み
//
// verified:
//   Yosupo Judge - Convolution (mod 1,000,000,007)
//     https://judge.yosupo.jp/problem/convolution_mod_1000000007
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// modint
//------------------------------//

// modint
template<int MOD = 998244353, bool PRIME = true> struct Fp {
    // inner value
    unsigned int val;
    
    // constructor
    constexpr Fp() : val(0) { }
    template<std::signed_integral T> constexpr Fp(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> constexpr Fp(T v) {
        val = (unsigned int)(v % get_umod());
    }
    constexpr long long get() const { return val; }
    constexpr static int get_mod() { return MOD; }
    constexpr static unsigned int get_umod() { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp() - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        return *this = *this * r.inv(); 
    }
    constexpr Fp pow(long long n) const {
        assert(n >= 0);
        Fp res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        assert(val);
        if (PRIME) {
            return pow(get_umod() - 2);
        } else {
            assert(gcd(val, get_umod()) == 1);
            long long m = get_umod(), a = val, b = m, u = 1, v = 0;
            while (b > 0) {
                auto t = a / b;
                a -= t * b, swap(a, b);
                u -= t * v, swap(u, v);
            }
            return Fp(u);
        }
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const Fp &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const Fp &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const Fp &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const Fp &r) const {
        return this->val >= r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp &x) {
        return os << x.val;
    }
    friend constexpr Fp pow(const Fp &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp inv(const Fp &r) {
        return r.inv();
    }
};


//------------------------------//
// NTT
//------------------------------//

// calc primitive root
constexpr int calc_primitive_root(long long m) {
    if (m == 1) return -1;
    if (m == 2) return 1;
    if (m == 998244353) return 3;
    if (m == 167772161) return 3;
    if (m == 469762049) return 3;
    if (m == 754974721) return 11;
    if (m == 645922817) return 3;
    if (m == 897581057) return 3;

    auto mod_pow = [&](long long a, long long n, long long m) {
        long long res = 1;
        while (n > 0) {
            if (n % 2 == 1) res = res * a % m;
            a = a * a % m;
            n >>= 1;
        }
        return res;
    };
    long long divs[20] = {};
    divs[0] = 2;
    long long cnt = 1;
    long long x = (m - 1) / 2;
    while (x % 2 == 0) x /= 2;
    for (long long i = 3; i * i <= x; i += 2) {
        if (x % i == 0) {
            divs[cnt++] = i;
            while (x % i == 0) x /= i;
        }
    }
    if (x > 1) divs[cnt++] = x;
    for (long long g = 2; ; g++) {
        bool ok = true;
        for (int i = 0; i < cnt; i++) {
            if (mod_pow(g, (m - 1) / divs[i], m) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
}

// NTT setup
template<class mint, int MOD = mint::get_mod(), int g = calc_primitive_root(mint::get_mod())>
struct ntt_setup {
    static constexpr int bsf_constexpr(unsigned int x) {
        int i = 0;
        while (!(x & (1 << i))) i++;
        return i;
    };

    static constexpr int rank = bsf_constexpr(MOD - 1);
    array<mint, rank + 1> root, iroot;  // root[i]^(2^i) = 1, root[i] * iroot[i] = 1
    array<mint, max(0, rank - 1)> rate2, irate2;
    array<mint, max(0, rank - 2)> rate3, irate3;

    ntt_setup() {
        root[rank] = mint(g).pow((MOD - 1) >> rank);
        iroot[rank] = root[rank].inv();
        for (int i = rank - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
        mint prod = 1, iprod = 1;
        for (int i = 0; i < rank - 1; i++) {
            rate2[i] = root[i + 2] * prod;
            irate2[i] = iroot[i + 2] * iprod;
            prod *= iroot[i + 2];
            iprod *= root[i + 2];
        }
        prod = 1, iprod = 1;
        for (int i = 0; i < rank - 2; i++) {
            rate3[i] = root[i + 3] * prod;
            irate3[i] = iroot[i + 3] * iprod;
            prod *= iroot[i + 3];
            iprod *= root[i + 3];
        }
    }
};

// NTT transformation
template<class mint, int MOD = mint::get_mod()> 
void ntt_trans(vector<mint> &v) {
    int n = (int)v.size();
    int h = 0;
    while ((1U << h) < (unsigned int)(n)) h++;
    static const ntt_setup<mint> setup;

    int len = 0;
    while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = v[i + offset];
                    auto r = v[i + offset + p] * rot;
                    v[i + offset] = l + r;
                    v[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len)) {
                    rot *= setup.rate2[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len++;
        } else {
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = setup.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot, rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * MOD * MOD;
                    auto a0 = 1ULL * v[i + offset].val;
                    auto a1 = 1ULL * v[i + offset + p].val * rot.val;
                    auto a2 = 1ULL * v[i + offset + p * 2].val * rot2.val;
                    auto a3 = 1ULL * v[i + offset + p * 3].val * rot3.val;
                    auto tmp = 1ULL * mint(a1 + mod2 - a3).val * imag.val;
                    auto na2 = mod2 - a2;
                    v[i + offset] = a0 + a2 + a1 + a3;
                    v[i + offset + p] = a0 + a2 + (mod2 * 2 - (a1 + a3));
                    v[i + offset + p * 2] = a0 + na2 + tmp;
                    v[i + offset + p * 3] = a0 + na2 + (mod2 - tmp);
                }
                if (s + 1 != (1 << len)) {
                    rot *= setup.rate3[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len += 2;
        }
    }
}

// NTT inv-transformation
template<class mint, int MOD = mint::get_mod()> 
void ntt_trans_inv(vector<mint> &v) {
    int n = (int)v.size();
    int h = 0;
    while ((1U << h) < (unsigned int)(n)) h++;
    static const ntt_setup<mint> setup;

    int len = h;
    while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = v[i + offset];
                    auto r = v[i + offset + p];
                    v[i + offset] = l + r;
                    v[i + offset + p] = (unsigned long long)((long long)(MOD) + l.val - r.val) * irot.val;
                }
                if (s + 1 != (1 << (len - 1))) {
                    irot *= setup.irate2[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len--;
        } else {
            int p = 1 << (h - len);
            mint irot = 1, iimag = setup.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot, irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * v[i + offset].val;
                    auto a1 = 1ULL * v[i + offset + p].val;
                    auto a2 = 1ULL * v[i + offset + p * 2].val;
                    auto a3 = 1ULL * v[i + offset + p * 3].val;
                    auto tmp = 1ULL * mint((MOD + a2 - a3) * iimag.val).val;
                    v[i + offset] = a0 + a1 + a2 + a3;
                    v[i + offset + p] = (a0 + (MOD - a1) + tmp) * irot.val;
                    v[i + offset + p * 2] = (a0 + a1 + (MOD - a2) + (MOD - a3)) * irot2.val;
                    v[i + offset + p * 3] = (a0 + (MOD - a1) + (MOD - tmp)) * irot3.val;
                }
                if (s + 1 != (1 << (len - 2))) {
                    irot *= setup.irate3[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len -= 2;
        }
    }
    mint in = mint(n).inv();
    for (int i = 0; i < n; i++) v[i] *= in;
}

// naive convolution
template<class VEC> VEC convolution_naive(const VEC &a, const VEC &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    VEC res(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; j++) for (int i = 0; i < n; i++) res[i + j] += a[i] * b[j];
    } else {
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) res[i + j] += a[i] * b[j];
    }
    return res;
}

// ntt convolution
template<class mint>
vector<mint> convolution_ntt(vector<mint> a, vector<mint> b) {
    int MOD = mint::get_mod();
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    assert((MOD - 1) % z == 0);
    a.resize(z), b.resize(z);
    ntt_trans(a), ntt_trans(b);
    for (int i = 0; i < z; i++) a[i] *= b[i];
    ntt_trans_inv(a);
    a.resize(n + m - 1);
    return a;
}

// convolution long long (if u64 is necessary, use convolution_ull)
template<class VEC> VEC convolution_ll(const VEC &a, const VEC &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return VEC();
    if (min(n, m) <= 60) return convolution_naive(a, b);

    static constexpr int MOD0 = 754974721;  // 2^24
    static constexpr int MOD1 = 167772161;  // 2^25
    static constexpr int MOD2 = 469762049;  // 2^26
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    vector<mint0> a0(n, 0), b0(m, 0);
    vector<mint1> a1(n, 0), b1(m, 0);
    vector<mint2> a2(n, 0), b2(m, 0);
    for (int i = 0; i < n; ++i) a0[i] = a[i], a1[i] = a[i], a2[i] = a[i];
    for (int i = 0; i < m; ++i) b0[i] = b[i], b1[i] = b[i], b2[i] = b[i];
    auto c0 = convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = convolution_ntt(std::move(a2), std::move(b2));

    VEC res(n + m - 1);
    long long mod0 = MOD0, mod01 = mod0 * MOD1;
    for (int i = 0; i < n + m - 1; ++i) {
        unsigned int y0 = c0[i].val;
        unsigned int y1 = (imod0 * (c1[i] - mint1(y0))).val;
        unsigned int y2 = (imod01 * (c2[i] - mint2(y0)) - imod1 * y1).val;
        res[i] = mod01 * y2 + mod0 * y1 + y0;
    }
    return res;
}

// convolution in general mod
template<class mint>
vector<mint> convolution_general_mod(const vector<mint> &a, const vector<mint> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return convolution_naive(a, b);
    if constexpr (std::is_same_v<mint, Fp<998244353>>) return convolution_ntt(a, b);

    static constexpr int MOD0 = 754974721;  // 2^24
    static constexpr int MOD1 = 167772161;  // 2^25
    static constexpr int MOD2 = 469762049;  // 2^26
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    vector<mint0> a0(n, 0), b0(m, 0);
    vector<mint1> a1(n, 0), b1(m, 0);
    vector<mint2> a2(n, 0), b2(m, 0);
    for (int i = 0; i < n; ++i) a0[i] = a[i].val, a1[i] = a[i].val, a2[i] = a[i].val;
    for (int i = 0; i < m; ++i) b0[i] = b[i].val, b1[i] = b[i].val, b2[i] = b[i].val;
    auto c0 = convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = convolution_ntt(std::move(a2), std::move(b2));

    vector<mint> res(n + m - 1);
    mint mod0 = MOD0, mod01 = mod0 * MOD1;
    for (int i = 0; i < n + m - 1; ++i) {
        unsigned int y0 = c0[i].val;
        unsigned int y1 = (imod0 * (c1[i] - mint1(y0))).val;
        unsigned int y2 = (imod01 * (c2[i] - mint2(y0)) - imod1 * y1).val;
        res[i] = mod01 * y2 + mod0 * y1 + y0;
    }
    return res;
}

// convolution overall
template<class T> vector<T> convolution(const vector<T> &a, const vector<T> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return convolution_naive(a, b);
    if constexpr (std::is_same_v<T, Fp<998244353>>) return convolution_ntt(a, b);
    else if constexpr (std::is_integral_v<T>) return convolution_ll(a, b);
    else return convolution_general_mod(a, b);
}


//------------------------------//
// Examples
//------------------------------//

// Yosupo Judge - Convolution (mod 1,000,000,007)
void Yosupo_Convolution_1000000007() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    
    int N, M;
    cin >> N >> M;

    vector<mint> a(N), b(M);
    for (int i = 0; i < N; ++i) cin >> a[i];
    for (int i = 0; i < M; ++i) cin >> b[i];
    
    auto c = convolution(a, b);
    
    for (int i = 0; i < c.size(); ++i) cout << c[i] << " ";
    cout << endl;
}


int main() {
    Yosupo_Convolution_1000000007();
}