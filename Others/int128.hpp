#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// int 128
using i128 = __int128;
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
        d--;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        d--;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
template <class T, enable_if_t<is_same_v<T, i128>, int> = 0>
constexpr i128 abs(T x) {
  return x < 0 ? -x : x;
}
constexpr i128 gcd(i128 a, i128 b) {
  while (b != 0) {
    i128 c = a % b;
    a = b, b = c;
  }
  return abs(a);
}

// uint 128
using u128 = __uint128_t;
u128 to_uinteger(const string &s) {
    u128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    return res;
}
istream& operator >> (istream &is, u128 &x) {
    string s;
    is >> s;
    x = to_uinteger(s);
    return is;
}
ostream& operator << (ostream &os, const u128 &x) {
    u128 ax = x;
    char buffer[128];
    char *d = end(buffer);
    do {
        d--;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        d--;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
