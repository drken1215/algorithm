//
// 削除可能な両端 priority queue
//
// verified:
//   Yosupo Library Checker - Double-Ended Priority Queue
//     https://judge.yosupo.jp/problem/double_ended_priority_queue
//

/*
 以下の機能をサポートする
 　・要素 x を挿入する
 　・要素 x を削除する (x が含まれていないとバグるので注意)
 　・最小値を取得する
 　・最小値を pop する
 　・最大値を取得する
 　・最大値を pop する
*/


#include <bits/stdc++.h>
using namespace std;


// Double Ended Priority Queue
template<class T> struct DoubleEndedPriorityQueue {
    // removable heap
    template<class QUETYPE> struct RemovablePriorityQueue {
        using VALTYPE = typename QUETYPE::value_type;
        QUETYPE que, delay;
        
        // constructor
        RemovablePriorityQueue() {}
        RemovablePriorityQueue(const RemovablePriorityQueue&) = default;
        RemovablePriorityQueue& operator = (const RemovablePriorityQueue&) = default;
        
        // getter
        constexpr int size() const { return (int)que.size() - (int)delay.size(); }
        constexpr bool empty() const { return size() == 0; }

        // push(x), remove(x)
        constexpr void push(VALTYPE x) { que.push(x); }
        constexpr void remove(VALTYPE x) { delay.push(x); }
        
        // pop min/max value
        constexpr VALTYPE pop() {
            T res = get();
            que.pop();
            return res;
        }
        
        // get min/max value (not pop)
        constexpr VALTYPE get() {
            assert(!que.empty());
            while (!delay.empty() && que.top() == delay.top()) {
                que.pop();
                delay.pop();
            }
            assert(!que.empty());
            return que.top();
        }
    };
    
    // inner data
    RemovablePriorityQueue<priority_queue<T, vector<T>, greater<T>>> min_que;
    RemovablePriorityQueue<priority_queue<T>> max_que;
    
    // constructor
    DoubleEndedPriorityQueue() {}
    DoubleEndedPriorityQueue(const DoubleEndedPriorityQueue&) = default;
    DoubleEndedPriorityQueue& operator = (const DoubleEndedPriorityQueue&) = default;
    
    // getter
    constexpr int size() const {
        return (int)min_que.size();
    }
    constexpr bool empty() const {
        return size() == 0;
    }
    
    // push(x), remove(x)
    constexpr void push(T x) {
        min_que.push(x);
        max_que.push(x);
    }
    constexpr void remove(T x) {
        min_que.remove(x);
        max_que.remove(x);
    }
    
    // get min, pop min
    constexpr T get_min() {
        return min_que.get();
    }
    constexpr T pop_min() {
        T x = min_que.pop();
        max_que.remove(x);
        return x;
    }
    
    // get max, pop max
    constexpr T get_max() {
        return max_que.get();
    }
    constexpr T pop_max() {
        T x = max_que.pop();
        min_que.remove(x);
        return x;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Double-Ended Priority Queue
void YosupoDoubleEndedProrityQueue() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    
    int N, Q;
    cin >> N >> Q;
    DoubleEndedPriorityQueue<int> que;
    for (int i = 0; i < N; ++i) {
        int v;
        cin >> v;
        que.push(v);
    }
    
    while (Q--) {
        int t;
        cin >> t;
        if (t == 0) {
            // x の挿入
            int x;
            cin >> x;
            que.push(x);
        } else if (t == 1) {
            // 最小の要素を取り出して除去する
            cout << que.pop_min() << endl;
        } else {
            // 最大の要素を取り出して除去する
            cout << que.pop_max() << endl;
        }
    }
}


int main() {
    YosupoDoubleEndedProrityQueue();
}

