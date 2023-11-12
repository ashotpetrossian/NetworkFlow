/**
 * @file networkFlowSolver.h
 * @author Ashot Petrossian (ashotpetrossian91@gmail.com)
 * @brief TODO:
 * @version 0.1
 * @date 2023-06-21
 * 
 * @copyright Copyright (c) 2023
 * 
 * @example TODO: 
 *
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <queue>
#include <memory>
#include <cmath>

namespace mySTL {

class Edge {
public:
    Edge(int from, int to, long capacity);
    Edge(const Edge& rhs);
    Edge& operator=(const Edge& rhs);
    Edge(Edge&& rhs) noexcept;
    Edge& operator=(Edge&& rhs) noexcept;
    ~Edge() = default;

    bool isResidual() const noexcept;
    long remainingCapacity() const noexcept;
    void augmentPath(long bottleneck);

    void setResidualEdge(const std::shared_ptr<Edge>& edge) noexcept;

    const int getToNode() const noexcept { return m_to; }
    const int getFromNode() const noexcept { return m_from; }
    
    friend std::ostream& operator<<(std::ostream& os, const Edge& edge) {
        os << "from " << edge.m_flow << " to " << edge.m_to
           << ", flow = " << edge.m_flow << ", capacity = " << edge.m_capacity
           << ", is residual: " << std::boolalpha << edge.isResidual() << std::endl;

        return os;
    }

private:
    int m_from;
    int m_to;
    long m_flow;
    long m_capacity;
    std::shared_ptr<Edge> m_residual;
};

namespace MiceAndOwlsProblem {

class MiceAndOwls {
    struct Mouse {
        std::pair<int, int> m_point;
        
        Mouse(int x, int y) : m_point{x, y} { }
    };

    struct Hole {
        int m_capacity;
        std::pair<int, int> m_point;

        Hole(int x, int y, int cap) : m_capacity{cap}, m_point{x, y} { }
    };

public:
    /// @brief for simplicity define your example here
    void executeMiceAndHolesSolution();

    void solveMicesAndOwlsProblem(const std::vector<std::unique_ptr<Mouse>>& mice, const std::vector<std::unique_ptr<Hole>>& holes, int radius);

    inline double calculateDistance(int x1, int y1, int x2, int y2) const noexcept {
        return std::sqrt(std::pow((x2 - x1), 2) + std::pow((y2 - y1), 2));
    }
};

} // MiceAndOwlsProblem

/**
 * @brief Abstract base class for solving network flow problem
 * Extenders will represent separate algorithm techniques.
 * 
 */
class NetworkFlowSolver {
public:
    NetworkFlowSolver(int vertexNum, int source, int sink);

    virtual void addEdge(int from, int to, long capacity);

    const std::vector<std::vector<std::shared_ptr<Edge>>>& getGraph() const noexcept { return m_graph; }
    long getMaxFlow() const noexcept { return m_maxFlow; }
    std::vector<int>& getVisitedArray() noexcept { return m_visited; } // extenders can modify the array
    int getSource() const noexcept { return m_source; }
    int getSink() const noexcept { return m_sink; }
    int& getVisitedToken() { return m_visitedToken; } // extenders can modity the value
    int getVertexNum() const noexcept { return m_vertexNum; }

    void visit(int index) { m_visited[index] = m_visitedToken; }
    bool isVisited(int vertex) { return m_visited[vertex] == m_visitedToken; }
    void markAllNodesAsVisited() { ++m_visitedToken; }
    
    void increaseMaxFlow(long bottleneck) { m_maxFlow += bottleneck; }

    /**
     * @brief this method should be overriden by the classes-extenders
     *  
     */
    virtual void execute() = 0;
private:
    int m_vertexNum;
    int m_source;
    int m_sink;

    int m_visitedToken {1}; // to avoid resetting visited array to false
    std::vector<int> m_visited;
    long m_maxFlow {0};

    std::vector<std::vector<std::shared_ptr<Edge>>> m_graph;
};

/**
 * @brief Ford-Fulkerson method using DFS 
 * Time complexity: O(E * maxFlow)
 * 
 */
class FordFulkersonDFSSolver : public NetworkFlowSolver {
public:
    FordFulkersonDFSSolver(int vertexNum, int source, int sink) : NetworkFlowSolver(vertexNum, source, sink) { }

    virtual void execute() override;
    long dfs(int source, long flow);
};

/**
 * @brief Edmonds-Karp method using BFS
 * Time complexity: O(V * E ^ 2)
 * Advanatges: Finds augmenting paths with minimum lengths, thus preventing to get lower bottleneck.
 *             The time complexity is independant from maxFlow
 * 
 */

class EdmondsKarpBFSSolver : public NetworkFlowSolver {
public:
    EdmondsKarpBFSSolver(int vertexNum, int source, int sink) : NetworkFlowSolver(vertexNum, source, sink) { }
    ~EdmondsKarpBFSSolver() = default;

    virtual void execute() override;
    long bfs();
};

/**
 * @brief Capacity-Scaling method used to prioritize edges with larger
 * capacities to avoid ending up with a path that has a small bottleneck.
 * 
 * The Capacity Scaling heuristic says that we should only take edges
 * whose remaining capacity is > delta (the first power of 2 of the largest capacity in the graph)
 * delta's step is /= 2
 * 
 * Time complexity: O(E ^ 2 * (log(U))), U is the biggest capacity edge in the graph
 * 
 */

class CapacityScalingSolver : public NetworkFlowSolver {
public:
    CapacityScalingSolver(int vertexNum, int source, int sink) : NetworkFlowSolver(vertexNum, source, sink) { }
    ~CapacityScalingSolver() = default;

    void addEdge(int from, int to, long capacity) override;
    virtual void execute() override;
    long dfs(int source, long flow);
    long& getDelta() { return m_delta; }
private:
    long m_delta = 0;
};

/**
 * @brief Dinic's algorithm
 * 
 * Main idea: Guide augmenting paths from s -> t using a level graph
 * The way Dinic's determine which edges make the progress towards the sink
 * and which do not is by building level graph. For levels we do BFS.
 * 
 * Time complexity: O(V ^ 2 * E)
 * For bipartite graphs O (sqrt(V * E))
 * 
 */

class DinicsSolver : public NetworkFlowSolver {
public:
    DinicsSolver(int vertexNum, int source, int sink) : NetworkFlowSolver(vertexNum, source, sink) {
        m_levels.resize(vertexNum);
    }

    virtual void execute();
    bool bfs();
    long dfs(int source, std::vector<int>& next, long flow);
    std::vector<int>& getLevels() { return m_levels; }
private:
    std::vector<int> m_levels;
};

} // mySTL