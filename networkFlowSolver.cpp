#include "networkFlowSolver.h"

using namespace mySTL;

Edge::Edge (int from, int to, long capacity) : m_from{from},
                                               m_to{to},
                                               m_capacity{capacity}
{ }

Edge::Edge(const Edge& rhs) : m_from{rhs.m_from},
                              m_to{rhs.m_to},
                              m_flow{rhs.m_flow},
                              m_capacity{rhs.m_capacity},
                              m_residual{rhs.m_residual}
{ }

Edge& Edge::operator=(const Edge& rhs)
{
    if (this != &rhs) {
        m_from = rhs.m_from;
        m_to = rhs.m_to;
        m_flow = rhs.m_flow;
        m_capacity = rhs.m_capacity;
        m_residual = rhs.m_residual;
    }

    return *this;
}

Edge::Edge(Edge&& rhs) noexcept : m_from{std::move(rhs.m_from)},
                                  m_to{std::move(rhs.m_to)},
                                  m_flow{std::move(rhs.m_flow)},
                                  m_capacity{std::move(rhs.m_capacity)},
                                  m_residual{std::move(rhs.m_residual)}
{
    rhs.m_residual.reset();
}

Edge& Edge::operator=(Edge&& rhs) noexcept
{
    if (this != &rhs) {
        m_from = std::move(rhs.m_from);
        m_to = std::move(rhs.m_to);
        m_flow = std::move(rhs.m_flow);
        m_capacity = std::move(rhs.m_capacity);
        m_residual = std::move(rhs.m_residual);
        rhs.m_residual.reset();
    }

    return *this;
}

bool Edge::isResidual() const noexcept
{
    return m_capacity == 0;
}

long Edge::remainingCapacity() const noexcept
{
    return m_capacity - m_flow;
}

void Edge::augmentPath(long bottleneck)
{
    m_flow += bottleneck;
    m_residual->m_flow -= bottleneck;
}

void Edge::setResidualEdge(const std::shared_ptr<Edge>& edge) noexcept
{
    m_residual = edge;
}

/**
 * @brief Construct a new Network Flow Solver
 * 
 * @param vertexNum - number of vertices
 * @param source - flow source
 * @param sink - flow sink
 */
NetworkFlowSolver::NetworkFlowSolver(int vertexNum, int source, int sink) : m_vertexNum{vertexNum},
                                                                            m_source{source},
                                                                            m_sink{sink}
{
    m_graph.resize(m_vertexNum);
    m_visited.resize(m_vertexNum, 0);
}

void NetworkFlowSolver::addEdge(int from, int to, long capacity)
{
    if (capacity < 0) throw std::logic_error("Forward edge cannot be negative\n");

    std::shared_ptr<Edge> edge1 = std::make_shared<Edge>(from, to, capacity);
    std::shared_ptr<Edge> edge2 = std::make_shared<Edge>(to, from, 0); // residual edge

    edge1->setResidualEdge(edge2);
    edge2->setResidualEdge(edge1);
    m_graph[from].push_back(edge1);
    m_graph[to].push_back(edge2);
}

/// @brief for debugging, should be deleted
std::vector<std::vector<int>> paths;
std::vector<int> flows;
std::vector<int> path;


/**
 * @brief runs dfs repeatedly till no augmenting path is found
 * the bottlenecks from each dfs is augmented to the maxFlow value
 */
void FordFulkersonDFSSolver::execute()
{
    for (long bottleneck = dfs(getSource(), INT32_MAX); bottleneck != 0; bottleneck = dfs(getSource(), INT32_MAX)) {
        ++getVisitedToken();
        increaseMaxFlow(bottleneck);
        flows.push_back(bottleneck);
        path.clear();
    }
    for (int i = 0; i < paths.size(); ++i) {
        std::cout << "Bottleneck: " << flows[i] << ": ";
        for (int v : paths[i]) {
            std::cout << v << " -> ";
        }
        std::cout << std::endl;
    }
}

/**
 * @brief dfs always starts from the source vertex of the graph
 * 
 * @param node - source vertex
 * @param flow - default big number. It will be reduced during the algorithm 
 *               while finging the smallest flow on the path
 * @return long - returns the bottleneck of the current path
 */
long FordFulkersonDFSSolver::dfs(int node, long flow)
{
    // std::cout << "considering vertex: " << node << " ";
    path.push_back(node);
    if (node == getSink()) {
        paths.push_back(path);
        return flow; // if the dfs reached to the sink
    }
    getVisitedArray()[node] = getVisitedToken(); // setting the current node as visited
    for (const std::shared_ptr<Edge>& edge : getGraph()[node]) {
        // if the current edge's flow is not the max and "to" vertex is not yet visited
        if (edge->remainingCapacity() > 0 && getVisitedArray()[edge->getToNode()] != getVisitedToken()) {
            // bottleneck value is the minimum of all remaining capacity values of all edges on the path
            // std::cout << "dfs on -> " << edge->getToNode() << std::endl;
            long bottleneck = dfs(edge->getToNode(), std::min(flow, edge->remainingCapacity()));
            // std::cout << "-> " << edge->getToNode() << " ";
            path.pop_back();

            // if we reached from source -> sink and the bottleneck is > 0
            // augment the flow with the bottleneck value
            if (bottleneck > 0) {
                // after backtracking if the bottleneck is not 0 we can augment the flow
                edge->augmentPath(bottleneck);
                return bottleneck;
            }
        }
    }

    return 0;
}

void EdmondsKarpBFSSolver::execute()
{
    for (long bottleneck = bfs(); bottleneck != 0; bottleneck = bfs()) {
        increaseMaxFlow(bottleneck);
        markAllNodesAsVisited();
    }
}

long EdmondsKarpBFSSolver::bfs()
{
    std::queue<int> q;
    int source = getSource();
    visit(source);
    q.push(source);

    std::vector<std::shared_ptr<Edge>> parents(getVertexNum(), nullptr);
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == getSink()) break;

        for (const std::shared_ptr<Edge>& edge : getGraph()[u]) {
            int v = edge->getToNode();
            if (edge->remainingCapacity() > 0 && !isVisited(v)) {
                visit(v);
                parents[v] = edge;

                q.push(v);
            }
        }
    }

    if (parents[getSink()] == nullptr) return 0; // the sink is no more reachable

    long bottleneck = INT32_MAX;

    /**
     * @brief iterating from sink to source by mapped incoming edges
     * 
     * @param edge != nullptr means that we reached the source as the source's parent edge is nullptr
     */
    std::cout << "Aumenting path: " << getSink() << " -> ";
    for (std::shared_ptr<Edge> edge = parents[getSink()]; edge != nullptr; edge = parents[edge->getFromNode()]) {
        std::cout << edge->getFromNode() << " -> ";
        bottleneck = std::min(bottleneck, edge->remainingCapacity());
    }
    std::cout << "\t\t bottleneck: " << bottleneck << std::endl;

    for (std::shared_ptr<Edge> edge = parents[getSink()]; edge != nullptr; edge = parents[edge->getFromNode()]) {
        edge->augmentPath(bottleneck);
    }

    return bottleneck;
}

void CapacityScalingSolver::addEdge(int from, int to, long capacity)
{
    NetworkFlowSolver::addEdge(from, to, capacity);
    getDelta() = std::max(getDelta(), capacity);
}

void CapacityScalingSolver::execute()
{
    getDelta() = static_cast<long>(std::pow(2, static_cast<long>(std::log2(getDelta()))));

    for (long bottleneck = 0; getDelta() > 0; getDelta() /= 2) {
        std::cout << "delta = " << getDelta() << std::endl;
        do {
            bottleneck = dfs(getSource(), INT32_MAX);
            increaseMaxFlow(bottleneck);
            std::cout << "bottleneck = " << bottleneck << std::endl;
            markAllNodesAsVisited();
            if (bottleneck > 0) flows.push_back(bottleneck);
            path.clear();
        } while (bottleneck != 0);
    }

    for (int i = 0; i < paths.size(); ++i) {
        std::cout << "Bottleneck: " << flows[i] << ": ";
        for (int v : paths[i]) {
            std::cout << v << " -> ";
        }
        std::cout << std::endl;
    }
}

long CapacityScalingSolver::dfs(int source, long flow)
{
    std::cout << "considering vertex: " << source << " -> ";
    path.push_back(source);
    if (source == getSink()) {
        paths.push_back(path);
        return flow;
    }
    visit(source);

    for (const std::shared_ptr<Edge>& edge : getGraph()[source]) {
        if (edge->remainingCapacity() > getDelta() && !isVisited(edge->getToNode())) {
            std::cout << "found valid toNode : " << edge->getToNode() << "  !!!" << std::endl;
            long bottleneck = dfs(edge->getToNode(), std::min(flow, edge->remainingCapacity()));
            
            path.pop_back();

            if (bottleneck > 0) {
                edge->augmentPath(bottleneck);
                return bottleneck;
            }
        }
    }
    return 0;
}

void DinicsSolver::execute()
{
    /// next[i] indicated the next edge index to take in the adjList for node i.
    /// This is part of Shimon Even and Alon Itai optimization for pruning dead ends as part of the DFS edges.
    std::vector<int> next(getVertexNum());

    while (bfs()) {
        for (long bottleneck = dfs(getSource(), next, INT32_MAX); bottleneck != 0; bottleneck = dfs(getSource(), next, INT32_MAX)) {
            increaseMaxFlow(bottleneck);
            flows.push_back(bottleneck);
            path.clear();
        }

        std::fill(next.begin(), next.end(), 0);
    }

    for (int i = 0; i < paths.size(); ++i) {
        std::cout << "path with bottleneck: " << flows[i] << " : ";
        for (int v : paths[i]) std::cout << v << " ";
        std::cout << std::endl;
    }
}

/// @brief levelize using bfs
/// @return if we were able to reach to the sink
bool DinicsSolver::bfs()
{
    std::vector<int>& levels {getLevels()};
    std::fill(levels.begin(), levels.end(), -1);

    std::queue<int> q;
    q.push(getSource());
    levels[getSource()] = 0;

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        for (const std::shared_ptr<Edge>& edge : getGraph()[node]) {
            if (edge->remainingCapacity() && levels[edge->getToNode()] == -1) {
                levels[edge->getToNode()] = levels[node] + 1; // mark levels
                q.push(edge->getToNode());
            }
        }
    }
    return levels[getSink()] != -1;
}

long DinicsSolver::dfs(int source, std::vector<int>& next, long flow)
{
    path.push_back(source);
    if (source == getSink()) {
        paths.push_back(path);
        return flow;
    }

    for (; next[source] < getGraph()[source].size(); ++next[source]) {
        const std::shared_ptr<Edge> edge = getGraph()[source][next[source]];
        if (edge->remainingCapacity() > 0 && getLevels()[edge->getToNode()] == getLevels()[source] + 1) {
            std::cout << "considering edge from: " << source << " to -> " << getGraph()[source][next[source]]->getToNode() << std::endl;
            long bottleneck = dfs(edge->getToNode(), next, std::min(flow, edge->remainingCapacity()));
            path.pop_back();

            if (bottleneck > 0) {
                edge->augmentPath(bottleneck);
                return bottleneck;
            }
        }
    }
    return 0;
}

void MiceAndOwlsProblem::MiceAndOwls::executeMiceAndHolesSolution()
{
    std::vector<std::unique_ptr<Mouse>> mice;
    mice.emplace_back(std::make_unique<Mouse>(1, 0));
    mice.emplace_back(std::make_unique<Mouse>(0, 1));
    mice.emplace_back(std::make_unique<Mouse>(8, 1));
    mice.emplace_back(std::make_unique<Mouse>(12, 0));
    mice.emplace_back(std::make_unique<Mouse>(12, 4));
    mice.emplace_back(std::make_unique<Mouse>(15, 5));

    std::vector<std::unique_ptr<Hole>> holes;
    holes.emplace_back(std::make_unique<Hole>(1, 1, 1));
    holes.emplace_back(std::make_unique<Hole>(10, 2, 2));
    holes.emplace_back(std::make_unique<Hole>(14, 5, 1));

    solveMicesAndOwlsProblem(mice, holes, /* radius = */ 3);   
}

void MiceAndOwlsProblem::MiceAndOwls::solveMicesAndOwlsProblem(const std::vector<std::unique_ptr<Mouse>>& mice,
                                                               const std::vector<std::unique_ptr<Hole>>& holes,
                                                               int radius)
{
    int m = mice.size();
    int h = holes.size();

    int n = m + h + 2;
    int s = n - 2;
    int t = n - 1;

    FordFulkersonDFSSolver solver(n, s, t);
    /// @brief connecting source to each mouse with capacity of one
    for (int i = 0; i < m; ++i) {
        solver.addEdge(s, i, 1);
    }

    /// @brief defining the edges mouse to hall (if the radius allows) with cap of 1
    for (int i = 0; i < m; ++i) {
        auto [mx, my] = mice[i]->m_point;
        for (int j = 0; j < h; ++j) {
            auto [hx, hy] = holes[j]->m_point;
            if (calculateDistance(mx, my, hx, hy) <= radius) {
                solver.addEdge(i, m + j, 1);
            }
        }
    }

    /// @brief defining the edges from holes to sink with cap of holes cap
    for (int i = 0; i < h; ++i) {
        solver.addEdge(m + i, t, holes[i]->m_capacity);
    }

    solver.execute();
    std::cout << "The max mice can be safe: " << solver.getMaxFlow() << std::endl;
}