#include <iostream>
#include <string>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigne-3.4.0/Eigen/Sparse>
#include <Snap-6.0/snap-core/Snap.h>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/global/ClusteringCoefficient.hpp>

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
    return A.sum() / (3 + 3 * static_cast<int>(!is_directed));
}

auto triangle_count_burkhardt_csr(int n, int m, bool is_directed) -> int {
    Eigen::SparseMatrix<int> A(n, n);

    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        A(u, v) += 1;
        A(v, u) += static_cast<int>(!is_directed);
    }

    A = (A * A).array() * A.array();
    return A.sum() / (3 + 3 * static_cast<int>(!is_directed));
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

auto triangle_count_sandia_сsr(int n, int m) -> int {
    Eigen::SparseMatrix A(n + 1, n + 1);

    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        if (u > v) std::swap(u, v);
        A(v, u) += 1;
    }

    A = (A * A).array() * A.array();
    return A.sum();
}


auto triangle_count_snap(int n, int m, bool is_directed) -> int {
    PUNGraph Graph = TUNGraph::New();
    for (int i = 0; i <= n; ++i) {
        Graph->AddNode(i);
    }
    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        Graph->AddEdge(u, v);
        if (is_directed) Graph->AddEdge(v, u);
    }

    return TSnap::GetTriads(Graph);
}

auto triangle_count_networkit(int n, int m, bool is_directed) -> int {
    NetworKit::Graph graph(n + 1);

    for (int _ = 0; _ < m; ++_) {
        int u, v;
        std::cin >> u >> v;
        graph.addEdge(u, v);
        graph.addEdge(v, u);
    }

    NetworKit::ClusteringCoefficient clusteringCoeff(graph);
    clusteringCoeff.run();

    return clusteringCoeff.numberOfTriangles();
}

auto main(int argc, char **argv) -> int {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <directed_flag> <vertices> <edges>\n";
        exit(1);
    }

    std::string algorithmName(argv[1]);
    bool isDirected = (std::string(argv[2]) == "1");
    int vertices = std::atoi(argv[3]);
    int edges = std::atoi(argv[4]);

    if (vertices <= 0 || edges < 0) {
        std::cerr << "Error: Number of vertices and edges must be positive.\n";
        exit(1);
    }

    if (algorithmName == "snap") {
        triangle_count_snap(vertices, edges, isDirected);
    } else if (algorithmName == "networkit") {
        triangle_count_networkit(vertices, edges, isDirected);
    } else if (algorithmName == "sandia") {
        triangle_count_sandia(vertices, edges);
    } else if (algorithmName == "sandia_csr") {
        triangle_count_sandia_сsr(vertices, edges);
    } else if (algorithmName == "burkhardt") {
        triangle_count_burkhardt(vertices, edges, isDirected);
    } else if (algorithmName == "burkhardt_csr") {
        triangle_count_burkhardt_csr(vertices, edges, isDirected);
    }

    std::cerr << "Error: Unknown algorithm '" << algorithmName << "'.\n";
    return 0;
}
