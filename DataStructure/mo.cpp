//
// Mo のアルゴリズム
//
// cf.
//   ei1333: Mo's algorithm
//     https://ei1333.hateblo.jp/entry/2017/09/11/211011
//
//
// verified (suffix array の lcp を sparse table で求める):
//   Yandex.Algorithm 2011 Round 2 D - Powerful array
//     http://codeforces.com/contest/86/problem/D
//


/*
 平方分割に基づいた、区間クエリに対する一般的なテク
 ・クエリ先読みができる、区間更新はなし
 ・区間の左端や右端を 1 個ずたしたものに対する値を高速に求められる (区間に a[idx] の要素を加えたり除いたり)
 
 ここでは、区間に含まれる数列の種類数を求める
 */


#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace std;

struct Mo {
    vector<int> left, right, index; // the interval's left, right, index
    vector<bool> v;
    int window;
    int nl, nr, ptr;
    
    Mo(int n) : window((int)sqrt(n)), nl(0), nr(0), ptr(0), v(n, false) { }
    
    /* push */
    void push(int l, int r) { left.push_back(l), right.push_back(r); }
    
    /* sort intervals */
    void build() {
        index.resize(left.size());
        iota(index.begin(), index.end(), 0);
        sort(begin(index), end(index), [&](int a, int b)
             {
                 if (left[a] / window != left[b] / window) return left[a] < left[b];
                 return bool((right[a] < right[b]) ^ (left[a] / window % 2));
             });
        
        //        sort(index.begin(), index.end(), [&](int a, int b) {
        //            if (left[a] / window != right[b] / window) return left[a] < left[b];
        //            else return right[a] < right[b];
        //        });
    }
    
    /* extend-shorten */
    void extend_shorten(int id) {
        v[id].flip();
        if (v[id]) insert(id);
        else erase(id);
    }
    
    /* next id of interval */
    int next() {
        if (ptr == index.size()) return -1;
        int id = index[ptr];
        while (nl > left[id]) extend_shorten(--nl);
        while (nr < right[id]) extend_shorten(nr++);
        while (nl < left[id]) extend_shorten(nl++);
        while (nr > right[id]) extend_shorten(--nr);
        return index[ptr++];
    }
    
    /* insert, erase (to be set appropriately) */
    void insert(int id);
    void erase(int id);
};



int N, Q;
int A[300001];
int res[200001];
int cnt[1000001];
int num_kind = 0;

void Mo::insert(int id) {
    int val = A[id];
    if (cnt[val] == 0) ++num_kind;
    ++cnt[val];
}

void Mo::erase(int id) {
    int val = A[id];
    --cnt[val];
    if (cnt[val] == 0) --num_kind;
}

int main() {
    scanf("%d", &N);
    for(int i = 0; i < N; i++) scanf("%d", &A[i]);
    scanf("%d", &Q);
    Mo mo(N);
    for(int i = 0; i < Q; i++) {
        int l, r; scanf("%d %d", &l, &r);
        mo.push(--l, r);
    }
    mo.build();
    for(int i = 0; i < Q; i++) res[mo.next()] = num_kind;
    for(int i = 0; i < Q; i++) printf("%d\n", res[i]);
}
