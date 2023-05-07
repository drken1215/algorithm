//
// 重心分解
//
// verified:
//   2018 第4回ドワンゴからの挑戦状 予選 E - ニワンゴくんの家探し
//     https://atcoder.jp/contests/dwacon2018-prelims/tasks/dwacon2018_prelims_e
//


#include <bits/stdc++.h>
using namespace std;


// sizeSubtree[v] := v を根とする部分ツリーのサイズ (分割統治の毎ステップごとに再利用)
// isRemoved[v] := v が既に取り除かれたかどうか
// whoIsParent[v] := ツリーDP時に v の親が誰だったか

using Graph = vector<vector<int>>;
struct TreeCenteroid {
    // input
    Graph tree;
    
    // results
    vector<int> centroids;
    
    // intermediate results
    vector<int> sizeSubtree, isRemoved, whoIsParent;
    
    // constructor
    TreeCenteroid() { }
    TreeCenteroid(const Graph &tree_) { init(tree_); }
    
    // init
    void init(const Graph &tree_) {
        tree = tree_;
        centroids.clear();
        sizeSubtree.resize((int)tree.size());
        whoIsParent.resize((int)tree.size());
        isRemoved.assign((int)tree.size(), false);
        for (int i = 0; i < (int)tree.size(); ++i) isRemoved[i] = false;
    }
    
    // subroutine
    void sub_find_centroid(int v, int size, int p = -1) {
        sizeSubtree[v] = 1;
        whoIsParent[v] = p;
        bool isCentroid = true;
        for (auto ch : tree[v]) {
            if (ch == p) continue;
            if (isRemoved[ch]) continue;
            sub_find_centroid(ch, size, v);
            if (sizeSubtree[ch] > size / 2) isCentroid = false;
            sizeSubtree[v] += sizeSubtree[ch];
        }
        if (size - sizeSubtree[v] > size / 2) isCentroid = false;
        if (isCentroid) centroids.push_back(v);
    }
    
    // first: centroid, second: vectors of (adj-node, size of adj-tree)
    pair<int, vector<pair<int,int>>> find_centroid(int root, int size) {
        vector<pair<int, int>> subtrees;
        centroids.clear();
        sub_find_centroid(root, size);
        int center = centroids[0];
        //isRemoved[center] = true;
        for (auto ch : tree[center]) {
            if (isRemoved[ch]) continue;
            if (ch == whoIsParent[center]) {
                subtrees.push_back(make_pair(ch, size - sizeSubtree[center]));
            }
            else {
                subtrees.push_back(make_pair(ch, sizeSubtree[ch]));
            }
        }
        return make_pair(center, subtrees);
    }
};


////////////////////////////////////////
// solver
////////////////////////////////////////

bool cmp(pair<int,int> a, pair<int,int> b) {
    swap(a.first, a.second);
    swap(b.first, b.second);
    return a > b;
}
 
void dwacon2018_prelims_E() {
    // 入力
    int N, Q;
    cin >> N >> Q;
    Graph tree(N);
    for (int i = 0; i < N - 1; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        tree[a].push_back(b);
        tree[b].push_back(a);
    }
    
    // 重心分解しながらクエリ処理
    int curnode = 0, cursize = N, res;
    TreeCenteroid tc(tree);
    while (Q--) {
        pair<int, vector<pair<int, int>>> c = tc.find_centroid(curnode, cursize);
        int g = c.first;
        vector<pair<int,int>> chs = c.second;
        sort(chs.begin(), chs.end(), cmp);
        if (chs.size() == 1) {
            cout << "? " << g + 1 << " " << chs[0].first + 1 << endl;
            int ans;
            cin >> ans;
            res = ans;
            break;
            
        } else {
            cout << "? " << chs[0].first + 1 << " " << chs[1].first + 1 << endl;
            int ans;
            cin >> ans;
            if (ans == chs[0].first + 1) {
                tc.isRemoved[g] = true;
                curnode = chs[0].first;
                cursize = chs[0].second;
                
            } else if (ans == chs[1].first + 1) {
                tc.isRemoved[g] = true;
                curnode = chs[1].first;
                cursize = chs[1].second;
                
            } else {
                tc.isRemoved[chs[0].first] = true;
                tc.isRemoved[chs[1].first] = true;
                curnode = g;
                cursize = cursize - chs[0].second - chs[1].second;
                
            }
            
            if (cursize == 1) {
                res = curnode + 1;
                break;
            }
        }
    }
    cout << "! " << res << endl;
}

int main() {
    dwacon2018_prelims_E();
}


