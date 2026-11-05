#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

using ll = long long;
using i128 = __int128_t;
using u128 = __uint128_t;
using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tll = array<long long, 3>;
using fll = array<long long, 4>;
using vint = vector<int>;
using vll = vector<long long>;
using dint = deque<int>;
using dll = deque<long long>;
using vvint = vector<vector<int>>;
using vvll = vector<vector<long long>>;
using vpll = vector<pair<long long, long long>>;
template<class T> using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
template<class S, class T> inline auto maxll(S a, T b) { return max(ll(a), ll(b)); }
template<class S, class T> inline auto minll(S a, T b) { return min(ll(a), ll(b)); }
template<class T> auto max(const T &a) { return *max_element(a.begin(), a.end()); }
template<class T> auto min(const T &a) { return *min_element(a.begin(), a.end()); }
template<class T> auto argmax(const T &a) { return max_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto argmin(const T &a) { return min_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto accum(const vector<T> &a) { return accumulate(a.begin(), a.end(), T()); }
template<class T> auto accum(const deque<T> &a) { return accumulate(a.begin(), a.end(), T()); }

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PF push_front
#define PB push_back
#define MP make_pair
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()
