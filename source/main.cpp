#include <iostream>
#include <string>
#include <eigen-3.4.0/Eigen/Dense>

#include "lib.hpp"

auto triangle_count_burkhardt(int n, int m, bool is_directed) -> int {
    Eigen::MatrixXd A(n + 1, n + 1);

    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        A(u, v) += 1;
        A(v, u) += static_cast<int>(!is_directed);
    }

    A = (A * A).array() * A.array();
    return A.sum() / 6;
}

auto triangle_count_sandia(int n, int m) -> int {
    Eigen::MatrixXd A(n + 1, n + 1);

    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        if (u > v) std::swap(u, v);
        A(v, u) += 1;
    }

    A = (A * A).array() * A.array();
    return A.sum();
}

auto main(int argc, char** argv) -> int {
    int n, m;
    std::cin >> n >> m;
    std::cout << triangle_count_sandia(n, m);
}
