#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// timer
struct Timer {
    bool started = false;
    chrono::high_resolution_clock::time_point start_time;

    Timer(bool started_ = false) : started(started_) {
        if (started_) start();
    }
    void start() {
        started = true;
        start_time = chrono::high_resolution_clock::now();
    }
    void stop() {
        started = false;
    }

    // second from start
    double operator () () {
        assert(started);
        chrono::high_resolution_clock::time_point target_time = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = target_time - start_time;
        return diff.count();
    }
};
