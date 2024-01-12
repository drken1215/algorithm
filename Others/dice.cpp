//
// Dice
//
// verified:
//   AOJ 0502 Dice
//   AOJ 1057 Rolling Dice
//


#include <iostream>
#include <vector>
using namespace std;


const vector<int> dx = {1, 0, -1, 0};
const vector<int> dy = {0, 1, 0, -1};

#define roll_swap(x, a, b, c, d) swap(x.a, x.b), swap(x.b, x.c), swap(x.c, x.d);
struct Dice {
    int top, front, right, left, back, bottom;
    
    Dice(int to = 1, int fr = 2, int ri = 3, int le = 4, int ba = 5, int bo = 6) 
    : top(to), front(fr), right(ri), left(le), back(ba), bottom(bo) { }
    void roll_right() { roll_swap((*this), top, left, bottom, right); }
    void roll_left() { roll_swap((*this), top, right, bottom, left); }
    void roll_front() { roll_swap((*this), top, back, bottom, front); }
    void roll_back() { roll_swap((*this), top, front, bottom, back); }
    void roll_cw() { roll_swap((*this), back, left, front, right); }
    void roll_ccw() { roll_swap((*this), back, right, front, left); }
    void roll(int dir) {
        if (dir == 0) (*this).roll_front();
        if (dir == 1) (*this).roll_right();
        if (dir == 2) (*this).roll_back();
        if (dir == 3) (*this).roll_left();
    }
    friend bool operator < (const Dice &d1, const Dice &d2) {
        int vd1[] = {d1.top, d1.front, d1.right, d1.left, d1.back, d1.bottom};
        int vd2[] = {d2.top, d2.front, d2.right, d2.left, d2.back, d2.bottom};
        return vector<int>(vd1, vd1 + 6) < vector<int>(vd2, vd2 + 6);
    }
    friend ostream& operator << (ostream& s, const Dice& d) { 
        return s << '<' << d.top << ',' << d.front << ',' << d.right << '>'; 
    }
};

vector<Dice> all_rotate() {
    vector<Dice> res(24);
    Dice d;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; j++) {
            res[i * 4 + j] = d;
            d.roll_cw();
        }
        if (i & 1) d.roll_front();
        else d.roll_right();
    }
    return res;    
}



//------------------------------//
// Examples
//------------------------------//

int main() {

}
