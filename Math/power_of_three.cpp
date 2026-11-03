//
// 平衡三進法展開
//
// verified
//   TopCoder SRM 604 DIV1 Easy PowerOfThree
//     https://community.topcoder.com/stat?c=problem_statement&pm=12917
//


#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<int> PowerThree(long long n) {
    vector<int> res;
    while (n != 0) {
        int amari = (n % 3 + 3) % 3;
        if (amari == 1) res.push_back(1), --n;
        else if (amari == 2) res.push_back(-1), ++n;
        else res.push_back(0);
        n /= 3;
    }
    return res;
}


class PowerOfThree {
public:
    string ableToGet(long long x, long long y) {
        vector<int> vx = PowerThree(x);
        vector<int> vy = PowerThree(y);

        // 桁数を合わせる
        while (vx.size() < vy.size()) vx.push_back(0);
        while (vy.size() < vx.size()) vy.push_back(0);
        
        // 各桁比較
        for (int i = 0; i < vx.size(); ++i) {
            if (vx[i] == 0 && vy[i] == 0) return "Impossible";
            if (vx[i] != 0 && vy[i] != 0) return "Impossible";
        }
        return "Possible";
    }
};



//------------------------------//
// Examples
//------------------------------//

int main() {

}
