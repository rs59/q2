#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_map>

// Hash function for pairs (used for edgeWeights unordered_map)
    struct HashPair {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ h2;
        }
    };

class Graph {
public:
    // Constructor
    Graph() {}

    // Copy constructor
    Graph(const Graph& other) {
        vertices = other.vertices;
        adjacencyList = other.adjacencyList;
        edgeWeights = other.edgeWeights;
        coarserToFinerMappings = other.coarserToFinerMappings;
        edgesMappings = other.edgesMappings;
        verticesWeightsMapping = other.verticesWeightsMapping;
    }

    // Function to add a vertex to the graph with its weight
    void addVertex(int vertexID, double weight) {
        vertices[vertexID] = weight;
        adjacencyList[vertexID] = std::vector<int>();
    }

    // Function to add an edge between two vertices with its weight
    void addEdge(int source, int destination, double weight) {
        adjacencyList[source].push_back(destination);
        adjacencyList[destination].push_back(source);
        edgeWeights[{source, destination}] = weight;
        edgeWeights[{destination, source}] = weight;
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
    std::vector<int>& getNeighbors(int vertexID) {
        return adjacencyList[vertexID];
    }

    // Function to get the number of vertices in the graph
    int numVertices() const {
        return vertices.size();
    }

    // Function to get the number of edges in the graph
    int numEdges() const {
        return edgeWeights.size()/2;
    }

    int size() const {
        return vertices.size();
    }

    //Function to get the vertices list and weights
    std::unordered_map<int, double> getVertices() {
        return vertices;
    }

    //Function to get the edges list and weights
    std::unordered_map<std::pair<int, int>, double, HashPair> getEdgeWeights(){
        return edgeWeights;
    }


    int getCoarsingLevel(){
        return coarserToFinerMappings.size();
    }

    void pushBackMapping(std::unordered_map<int,int>mapping){
        coarserToFinerMappings.push_back(mapping);
    }

    std::unordered_map<int,int> getMapping(int level){
        return coarserToFinerMappings[level];
    }

    void pushBackMappingEdges(std::unordered_map<std::pair<int, int>, double, HashPair> edges){
        edgesMappings.push_back(edges);
    }

    std::unordered_map<std::pair<int, int>, double, HashPair> getMappingEdges(int level){
        return edgesMappings[level];
    }

    void pushBackVerticesWeights(std::unordered_map<int, double> weights){
        verticesWeightsMapping.push_back(weights);
    }

    std::unordered_map<int, double> getMappingVerticesWeights(int level){
        return verticesWeightsMapping[level];
    }

    void print(){
        std::cout << "Vertices: " << std::endl;
        for (const auto& entry : vertices) {
            int key = entry.first;
            double value = entry.second;
            std::cout << "Vertex: " << key << ", Weight: " << value << std::endl;
        }

        std::cout << "Edges: " << std::endl;
        for (const auto& entry : edgeWeights) {
            const auto& key = entry.first;
            double value = entry.second;
            std::cout << "Edge: (" << key.first << ", " << key.second << "), Weight: " << value << std::endl;
        }
    }

private:
    // Data structures to store the graph data
    std::unordered_map<int, double> vertices;  // VertexID -> Weight
    std::unordered_map<int, std::vector<int>> adjacencyList;  // VertexID -> List of neighbors
    std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights; // (Source, Destination) -> Edge Weight

    // Container to hold vertex and edges mappings for each level
    std::vector<std::unordered_map<int, int>> coarserToFinerMappings;
    std::vector<std::unordered_map<std::pair<int, int>, double, HashPair>> edgesMappings;
    std::vector<std::unordered_map<int, double>> verticesWeightsMapping;
};

#endif // GRAPH_H