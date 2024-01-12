//
// エラトスネテスの篩を用いた各種アルゴリズム
//　・高速素因数分解
//　・約数列挙
//　・メビウス関数
//
// cf.
//   高速素因数分解
//     http://www.osak.jp/diary/diary_201310.html#20131017
//
// verified
//   Codeforces 613 DIV2 F. Classical?
//     https://codeforces.com/contest/1285/problem/F
//


#include <iostream>
#include <vector>
#include <stack>
using namespace std;


// isprime[n] := is n prime?
// mebius[n] := mebius value of n
// min_factor[n] := the min prime-factor of n

struct Eratos {
    vector<int> primes;
    vector<bool> isprime;
    vector<int> mebius;
    vector<int> min_factor;

    Eratos(int MAX) : primes(),
                      isprime(MAX+1, true),
                      mebius(MAX+1, 1),
                      min_factor(MAX+1, -1) {
        isprime[0] = isprime[1] = false;
        min_factor[0] = 0, min_factor[1] = 1;
        for (int i = 2; i <= MAX; ++i) {
            if (!isprime[i]) continue;
            primes.push_back(i);
            mebius[i] = -1;
            min_factor[i] = i;
            for (int j = i*2; j <= MAX; j += i) {
                isprime[j] = false;
                if ((j / i) % i == 0) mebius[j] = 0;
                else mebius[j] = -mebius[j];
                if (min_factor[j] == -1) min_factor[j] = i;
            }
        }
    }

    // prime factorization
    vector<pair<int,int>> prime_factors(int n) {
        vector<pair<int,int> > res;
        while (n != 1) {
            int prime = min_factor[n];
            int exp = 0;
            while (min_factor[n] == prime) {
                ++exp;
                n /= prime;
            }
            res.push_back(make_pair(prime, exp));
        }
        return res;
    }

    // enumerate divisors
    vector<int> divisors(int n) {
        vector<int> res({1});
        auto pf = prime_factors(n);
        for (auto p : pf) {
            int n = (int)res.size();
            for (int i = 0; i < n; ++i) {
                int v = 1;
                for (int j = 0; j < p.second; ++j) {
                    v *= p.first;
                    res.push_back(res[i] * v);
                }
            }
        }
        return res;
    }
};



//------------------------------//
// Examples
//------------------------------//

long long GCD(long long x, long long y) {
    if (y == 0) return x;
    return GCD(y, x % y);
}

const int MAX = 110000;

int main() {
    Eratos er(MAX);
    vector<vector<int>> divs(MAX);
    for (int n = 1; n < MAX; ++n) divs[n] = er.divisors(n);

    stack<int> st;
    vector<int> counter(MAX, 0);
    auto push = [&](int x) {
        st.push(x);
        for (auto d : divs[x]) ++counter[d];
    };
    auto pop = [&]() {
        int x = st.top();
        st.pop();
        for (auto d : divs[x]) --counter[d];
    };
    auto count = [&](int x) {
        int res = 0;
        for (auto d : divs[x]) res += er.mebius[d] * counter[d];
        return res;
    };

    int N;
    cin >> N;
    vector<long long> a(N);
    vector<bool> isin(MAX+1, false);
    long long res = 0;
    for (int i = 0; i < N; ++i) {
        cin >> a[i];
        res = max(res, a[i]);
        isin[a[i]] = true;
    }

    for (int g = 1; g < MAX; ++g) {
        for (int v = (MAX/g)*g; v >= g; v -= g) {
            if (!isin[v]) continue;
            while (count(v/g)) {
                long long t = st.top() * g;
                res = max(res, t / GCD(v, t) * v);
                pop();
            }
            push(v/g);
        }
        while (!st.empty()) pop();
    }
    cout << res << endl;
}
