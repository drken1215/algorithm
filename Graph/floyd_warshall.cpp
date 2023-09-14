//
// 全頂点対間最短路を求める、Floyd-Warshall のアルゴリズム
//
// cf.
//
//
// verified:
//   AOJ Course GRL_1_C Shortest Path - All Pairs Shortest Path
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_C&lang=jp
//


/*
    dp[a][b] := a から　b への距離
 
 として初期化しておいて、アルゴリズム終了後には
 
    dp[a][b] := a から　b への最短距離
 
 が格納される。注意点として、初期化時に dp[a][a] = 0 となるようにする必要がある。
 なお負閉路判定に用いることもできて、
 
 - 負閉路がなければ、任意の頂点 a に対して dp[a][a] = 0
 - 負閉路があったら、ある頂点 a に対して dp[a][a] < 0
 
 */


#include <iostream>
#include <vector>
using namespace std;
const long long INF = 1LL<<60;

int main() {
    int V, E; cin >> V >> E;
    
    // initialization of Floyd-Warshall
    vector<vector<long long> > dp(V, vector<long long>(V, INF));
    for (int i = 0; i < V; ++i) dp[i][i] = 0;  // necessary
    
    // input
    for (int e = 0; e < E; ++e) {
        int a, b, w; cin >> a >> b >> w;
        dp[a][b] = w;
    }
    
    // Floyd-Warshall
    for (int k = 0; k < V; ++k)
        for (int i = 0; i < V; ++i)
            for (int j = 0; j < V; ++j)
                dp[i][j] = min(dp[i][j], dp[i][k] + dp[k][j]);
    
    // output
    bool isnegative = false;
    for (int v = 0; v < V; ++v) if (dp[v][v] < 0) isnegative = true;
    if (isnegative) puts("NEGATIVE CYCLE");
    else {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (dp[i][j] < INF/2) cout << dp[i][j];
                else cout << "INF";
                if (j != V-1) cout << " ";
            }
            cout << endl;
        }
    }
}
