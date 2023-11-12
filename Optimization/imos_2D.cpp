//
// imos 2D algorithm
//
// cf.
//   https://imoz.jp/algorithms/imos_method.html
//
// verified
//   AOJ Course DSL_5_B Cumulative Sum - The Maximum Number of Overlaps
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_5_B&lang=jp
//


#include <iostream>
#include <vector>
using namespace std;

const int MAX_X = 1001;
const int MAX_Y = 1001;

using Point = pair<int,int>;          // 2 ������ɸ
using Rectangle = pair<Point,Point>;  // Ĺ������ (����, ����)

int imos(const vector<Rectangle> &recs) {
    vector<vector<int> > nums(MAX_X+1, vector<int>(MAX_Y+1, 0));
    for (auto rec : recs) {
        nums[rec.first.first][rec.first.second]++;
        nums[rec.first.first][rec.second.second]--;
        nums[rec.second.first][rec.first.second]--;
        nums[rec.second.first][rec.second.second]++;
    }
    for (int x = 0; x < MAX_X; ++x) 
        for (int y = 0; y < MAX_Y; ++y)
            nums[x][y+1] += nums[x][y];
    for (int y = 0; y < MAX_Y; ++y)
        for (int x = 0; x < MAX_X; ++x)
            nums[x+1][y] += nums[x][y];
    int res = 0;
    for (int x = 0; x <= MAX_X; ++x)
        for (int y = 0; y <= MAX_Y; ++y)
            res = max(res, nums[x][y]);
    return res;
}


int main() {
    int N; cin >> N;
    vector<Rectangle> recs(N);
    for (int i = 0; i < N; ++i)
        cin >> recs[i].first.first >> recs[i].first.second >> recs[i].second.first >> recs[i].second.second;
    cout << imos(recs) << endl;
}
