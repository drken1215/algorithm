//
// 削除可能ヒープ
//
// reference:
//   AtCoder ABC 281 E - Least Elements の解説
//     https://drken1215.hatenablog.com/entry/2022/12/14/185200
//
// verified:
//   AtCoder ABC 281 E - Least Elements
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
//
//   AtCoder ABC 330 E - Mex and Update
//     https://atcoder.jp/contests/abc330/tasks/abc330_e
//


/*
 以下の処理ができる: 集合 S に対して
 　・要素 x を挿入する
 　・要素 x を削除する (x が S に含まれていないとバグる)
 　・最小値 (or 最大値) を取得する
*/


#include <bits/stdc++.h>
using namespace std;


// removable min heap
template<class T> struct removable_min_heap {
    // inner data
    priority_queue<T, vector<T>, greater<T>> que, delay;
    
    // constructor
    removable_min_heap() {}

    // add(x), remove(x)
    void add(T x) { que.push(x); }
    void remove(T x) { delay.push(x); }
    int size() { return (int)que.size() - (int)delay.size(); }
    
    // pop min value
    T pop() {
        T res = get_min();
        que.pop();
        return res;
    }
    
    // get min value (not pop)
    T get_min() {
        assert(!que.empty());
        while (!delay.empty() && que.top() == delay.top()) {
            que.pop();
            delay.pop();
        }
        assert(!que.empty());
        return que.top();
    }
};

// removable max heap
template<class T> struct removable_max_heap {
    // inner data
    priority_queue<T> que, delay;
    
    // constructor
    removable_max_heap() {}

    // add(x), remove(x), size
    void add(T x) { que.push(x); }
    void remove(T x) { delay.push(x); }
    int size() { return (int)que.size() - (int)delay.size(); }
    
    // pop min value
    T pop() {
        T res = get_max();
        que.pop();
        return res;
    }
    
    // get max value (not pop)
    T get_max() {
        assert(!que.empty());
        while (!delay.empty() && que.top() == delay.top()) {
            que.pop();
            delay.pop();
        }
        assert(!que.empty());
        return que.top();
    }
};



//------------------------------//
// Examples
//------------------------------//

// ABC 281 E - Least Elements
void ABC_281_E() {
    // 入力
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 小さい順に K 個の総和
    long long sum = 0;
    
    // 左側と右側の priority queue
    removable_max_heap<long long> left;
    removable_min_heap<long long> right;
    
    // add x
    auto add = [&](long long x) -> void {
        // とりあえず x を left に挿入する
        left.add(x);
        sum += x;

        // left のサイズが K を超えるなら left の最大値を right に移す
        if (left.size() > K) {
            long long y = left.pop();
            right.add(y);
            sum -= y;
        }
    };
    
    // remove x
    auto remove = [&](long long x) -> void {
        if (x <= left.get_max()) {
            // x が left の最大値以下ならば、x を left から削除する
            left.remove(x);
            sum -= x;
        } else {
            // x がそれより大きければ、x を right から削除する
            right.remove(x);
        }

        // left のサイズが K 未満になるなら right の最小値を left に移す
        if (left.size() < K) {
            long long y = right.pop();
            left.add(y);
            sum += y;
        }
    };
    for (int i = 0; i < M; ++i) add(A[i]);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << sum << " ";
        if (i+M < N) add(A[i+M]), remove(A[i]);
    }
    cout << endl;
}


// ABC 330 E - Mex and Update
void ABC_330_E() {
    // 入力
    int N, Q;
    cin >> N >> Q;
    vector<int> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    // データ
    map<int,int> ma;  // ma[val] := A の中に val が何個あるか
    removable_min_heap<int> que;
    for (int i = 0; i <= N; ++i) que.add(i);
    for (auto x : A) {
        if (ma[x] == 0) que.remove(x);
        ++ma[x];
    }
    
    // A の中から値 x を 1 個削除
    auto del = [&](int x) -> void {
        --ma[x];
        if (ma[x] == 0) que.add(x);
    };
    
    // A に値 x を 1 個挿入
    auto ins = [&](int x) -> void {
        if (ma[x] == 0) que.remove(x);
        ++ma[x];
    };
            
    // クエリ処理
    while (Q--) {
        int id, x;
        cin >> id >> x;
        --id;
        
        del(A[id]);
        A[id] = x;
        ins(A[id]);
        
        cout << que.get_min() << endl;
    }
}


int main() {
    //ABC_281_E();
    ABC_330_E();
}


