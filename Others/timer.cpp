#include <iostream>
#include <sys/time.h>
using namespace std;


double getTime() {
    struct timeval s;
    gettimeofday(&s, NULL);
    return s.tv_sec + s.tv_usec * 1e-6;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    auto s = getTime();
    for (int i = 0; i < 1000000000; ++i) {

    }
    auto t = getTime();
    cout << t - s << endl;
}
