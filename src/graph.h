#include <iostream>
#include <vector>
#include <unordered_map>

class Graph {
public:
    // Constructor
    Graph() {}

    // Function to add a vertex to the graph with its weight
    void addVertex(int vertexID, double weight) {
        vertices[vertexID] = weight;
        adjacencyList[vertexID] = std::vector<int>();
    }

    // Function to add an edge between two vertices with its weight
    void addEdge(int source, int destination, double weight) {
        adjacencyList[source].push_back(destination);
        edgeWeights[{source, destination}] = weight;
    }

    // Function to get the weight of a vertex
    double getVertexWeight(int vertexID) {
        return vertices[vertexID];
    }

    // Function to get the weight of an edge
    double getEdgeWeight(int source, int destination) {
        return edgeWeights[{source, destination}];
    }

    // Function to check if a vertex exists in the graph
    bool containsVertex(int vertexID) {
        return vertices.find(vertexID) != vertices.end();
    }

    // Function to check if an edge exists in the graph
    bool containsEdge(int source, int destination) {
        return edgeWeights.find({source, destination}) != edgeWeights.end();
    }

    // Function to get the neighbors of a vertex
    const std::vector<int>& getNeighbors(int vertexID) {
        return adjacencyList[vertexID];
    }

    // Function to get the number of vertices in the graph
    int numVertices() const {
        return vertices.size();
    }

    // Function to get the number of edges in the graph
    int numEdges() const {
        return edgeWeights.size();
    }

    //Function to get the vertices list and weights
    std::unordered_map<int, double> getVertices(){
        return vertices;
    }

    //Function to get the edges list and weights
    std::unordered_map<int, double> getEdgeWeights(){
        return edge_weights;
    }

private:
    // Data structures to store the graph data
    std::unordered_map<int, double> vertices;  // VertexID -> Weight
    std::unordered_map<int, std::vector<int>> adjacencyList;  // VertexID -> List of neighbors
    std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights; // (Source, Destination) -> Edge Weight

    // Hash function for pairs (used for edgeWeights unordered_map)
    struct HashPair {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ h2;
        }
    };
};