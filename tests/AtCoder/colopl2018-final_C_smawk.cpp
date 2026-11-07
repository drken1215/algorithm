#define PROBLEM "https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c"

#include "smawk.hpp"

int main() {
    long long N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    auto func = [&](long long i, long long j) -> long long {
        return A[j] + (j - i) * (j - i);
    };
    auto res = SMAWK<long long>(N, N, func);
    for (int i = 0; i < N; i++) cout << res[i].first << endl;
}
