//
// Mo のアルゴリズム
//
// cf.
//   ei1333: Mo's algorithm
//     https://ei1333.hateblo.jp/entry/2017/09/11/211011
//
//
// verified (suffix array の lcp を sparse table で求める):
//   SPOJ FREQUENT - Frequent values
//     https://www.spoj.com/problems/FREQUENT/
//


/*
 平方分割に基づいた、区間クエリに対する一般的なテク
 ・クエリ先読みができる、区間更新はなし
 ・区間の左端や右端を 1 個ずたしたものに対する値を高速に求められる (区間に a[idx] の要素を加えたり除いたり)
 
 ここでは
 cnt[i] := i が何個あるか
 hist[c] := 区間内に c 個ある値が何種類あるか
 num_kind := 区間内に何種類の数があるか
 mode := 最頻値の出現回数
 */


#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstring>
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


const int GETA = 100000;
int N, Q;
int A[100100], res[100100];
int cnt[200100], hist[100100];
int num_kind = 0, mode = 0;

void Mo::insert(int id) {
    int val = A[id];
    if (cnt[val] == 0) ++num_kind;
    --hist[cnt[val]];
    ++cnt[val];
    ++hist[cnt[val]];
    mode = max(mode, cnt[val]);
}

void Mo::erase(int id) {
    int val = A[id];
    --hist[cnt[val]];
    if(cnt[val] == mode && hist[cnt[val]] == 0) --mode;
    --cnt[val];
    ++hist[cnt[val]];
    if (cnt[val] == 0) --num_kind;
}

int main() {
    while (scanf("%d", &N), N) {
        memset(cnt, 0, sizeof(cnt));
        memset(hist, 0, sizeof(hist));
        mode = 0;
        
        scanf("%d", &Q);
        for(int i = 0; i < N; i++) {
            scanf("%d", &A[i]);
            A[i] += GETA;
        }
        Mo mo(N);
        for(int i = 0; i < Q; i++) {
            int x, y;
            scanf("%d %d", &x, &y);
            mo.push(--x, y);
        }
        mo.build();
        for(int i = 0; i < Q; i++) res[mo.next()] = mode;
        for(int i = 0; i < Q; i++) printf("%d\n", res[i]);
    }
}
