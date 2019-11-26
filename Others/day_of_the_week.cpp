//
// pn + r (n は非負整数) で表せる整数のうち、x 以上となる最小の整数
//
// verified:
//   Codeforces 304 DIV2 183 - B
//


#include <iostream>
using namespace std;

// グレゴリオ歴を求める
// グレゴリオ歴におけるユリウス日（1582年10月15日を基準とした日、但し、-4712年1月1日が0となるようにする）
int GregorianDay(int y, int m, int d) {
    y += 4800;
    if (m < 3) {--y; m += 12;}
    return 365*y+y/4-y/100+y/400+(153*m-457)/5+d-32083;
}

// ユリウス歴におけるユリウス日（-4712年1月1日を0日目とした日）
int JulianDay(int y, int m, int d) {
    y += 4716;
    if (m < 3) {--y; m += 12;}
    return 365*y+y/4+(153*m-457)/5+d-1402;
}

// グレゴリオ歴の下での曜日（日:1、月:2、火:3、水:4、木:5、金:6、土:0）
int dayWeek(int y, int m, int d) {
    return (GregorianDay(y, m, d) + 5) % 7;
}

int main() {
    cout << dayWeek(2019, 11, 27) << endl;
}
