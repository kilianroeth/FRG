#pragma once

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

// returns a linspace like numpy
template <typename T>
std::vector<T> linspace(T min, T max, size_t N) {
    if (max < min) {
        std::swap(min, max);
    }
    std::vector<T> linspace(N);
    if (N == 0) return linspace;
    for (size_t i = 0; i < N; ++i) {
        linspace[i] = min + i*(max - min)/(N - 1);
    }
    return linspace;
}

// return a logspace, from 10^min to 10^max with N10 values per step size
template <typename T>
std::vector<T> logspace(T min, T max, size_t N10) {
    if (max < min) {
        std::swap(min, max);
    }
    size_t N = (max - min) * N10;
    std::vector<T> logspace(N);
    if (N == 0) return logspace;
    for (size_t i = 0; i < N; ++i) {
        T exponent = min + static_cast<T>(i)*(max - min)/(N - 1);
        logspace[i] = std::pow(10, exponent);
    }
    return logspace;
}

// prints a progress bar
inline void progressBar(size_t current, size_t total, size_t width = 50, std::ostream& os = std::cerr) {
    if (total == 0) {
        os << "\r[" << std::string(width, '-') << "] 100%" << std::flush;
        return;
    }

    if (current > total) {
        current = total;
    }

    const size_t filled = (current * width) / total;
    const size_t percent = (current * 100) / total;

    os << '\r' << '[';
    for (size_t i = 0; i < width; ++i) {
        os << (i < filled ? '#' : '-');
    }
    os << "] " << percent << '%' << std::flush;

    if (current == total) {
        os << '\n';
    }
}

// compute the volume of a d sphere
inline double Ω(size_t d) {
    return 2.*pow(M_PI,d/2.) / std::tgamma(d/2.);
}