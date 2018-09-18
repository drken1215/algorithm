#include <iostream>
#include <string>
#include <vector>
using namespace std;

template<class T> struct CHT {
    struct Line {
        T a, b;
        Line(T a = 0, T b = 0) : a(a), b(b) { }
        T get(T x) {
            return (x - a) * (x - a) + b;
        }
    };
    
    struct Node {
        Line line;
        Node *left, *right;
        Node(Line line) : line(line), left(nullptr), right(nullptr){}
    };
    
    const T MIN, MAX, INF;
    Node* root;
    
    CHT(T MIN, T MAX, T INF) : MIN(MIN), MAX(MAX), INF(INF), root(nullptr) { }
    Node* insert(Node* p, T low, T high, Line& l){
        if(!p) return new Node(l);
        if (p->line.get(low) <= l.get(low) && p->line.get(high) <= l.get(high)) return p;
        if (p->line.get(low) >= l.get(low) && p->line.get(high) >= l.get(high)) {
            p->line = l;
            return p;
        }
        T mid = (low + high) / 2;
        if (p->line.get(mid) > l.get(mid)) swap(p->line , l);
        if (p->line.get(low) >= l.get(low)) p->left = insert(p->left, low, mid, l);
        else p->right = insert(p->right, mid, high, l);
        return p;
        
    }
    void insert(T a, T b){
        Line l(a, b);
        root = insert(root, MIN, MAX, l);
    }
    T query(Node* p, T low, T high, T x){
        if(!p) return INF;
        if(low == high) return p->line.get(x);
        T mid = (low + high) / 2;
        if (x <= mid) return min(p->line.get(x), query(p->left, low, mid, x));
        else return min(p->line.get(x), query(p->right, mid, high, x));
    }
    T query(T x){
        return query(root, MIN, MAX, x);
    }
};


const long long INF = 1LL<<60;
int main() {
    int N, K; cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    vector<vector<long long> > dp(N+1, vector<long long>(K+1, INF));
    dp[0][0] = 0;
    for (int k = 0; k < K; ++k) {
        CHT<long long> cht(0, 1100000, INF);
        for (int i = 0; i < N; ++i) {
            cht.insert(A[i], dp[i][k]);
            dp[i+1][k+1] = cht.query(A[i]);
        }
    }
    cout << dp[N][K] << endl;
}
