#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#define INF_EDGE_WEIGHT 2000000000.0

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
        this->clear();

        vertices = other.vertices;
        adjacencyList = other.adjacencyList;
        edgeWeights = other.edgeWeights;

        for(std::unordered_map<int, int> coarser_to_finer_mapping : other.coarserToFinerMappings){
            this->pushBackMapping(coarser_to_finer_mapping);
        }
        for(std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping : other.edgesMappings){
            this->pushBackMappingEdges(edgesMapping);
        }
        for(std::unordered_map<int, double> verticesWeights : other.verticesWeightsMapping){
            this->pushBackVerticesWeights(verticesWeights);
        }
    }

    // Copy assignment operator
    Graph& operator=(const Graph& other) {
        if (this == &other) {
            // Handle self-assignment
            return *this;
        }

        this->clear();

        vertices = other.vertices;
        adjacencyList = other.adjacencyList;
        edgeWeights = other.edgeWeights;

        for(std::unordered_map<int, int> coarser_to_finer_mapping : other.coarserToFinerMappings){
            this->pushBackMapping(coarser_to_finer_mapping);
        }
        for(std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping : other.edgesMappings){
            this->pushBackMappingEdges(edgesMapping);
        }
        for(std::unordered_map<int, double> verticesWeights : other.verticesWeightsMapping){
            this->pushBackVerticesWeights(verticesWeights);
        }

        return *this;
    }

    // Function to add a vertex to the graph with its weight
    void addVertex(int vertexID, double weight = 1) {
        vertices[vertexID] = weight;
        adjacencyList[vertexID] = std::vector<int>();
    }

    void copyCoarseningData(Graph& other){
        // Clear mappings and other data
        coarserToFinerMappings.clear();
        edgesMappings.clear();
        verticesWeightsMapping.clear();

        for(std::unordered_map<int, int> coarser_to_finer_mapping : other.coarserToFinerMappings){
            this->pushBackMapping(coarser_to_finer_mapping);
        }
        for(std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping : other.edgesMappings){
            this->pushBackMappingEdges(edgesMapping);
        }
        for(std::unordered_map<int, double> verticesWeights : other.verticesWeightsMapping){
            this->pushBackVerticesWeights(verticesWeights);
        }
    }

    // Function to add an edge between two vertices with its weight
    void addEdge(int source, int destination, double weight) {
        // Check if the destination is not already in the adjacency list of the source
        if (std::find(adjacencyList[source].begin(), adjacencyList[source].end(), destination) == adjacencyList[source].end()) {
            adjacencyList[source].push_back(destination);
        }

        // Check if the source is not already in the adjacency list of the destination
        if (std::find(adjacencyList[destination].begin(), adjacencyList[destination].end(), source) == adjacencyList[destination].end()) {
            adjacencyList[destination].push_back(source);
        }

        // Set the edge weight
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

        // Function to get the weight of an edge
    double getEdgeWeightKL(const int& source, const int& destination) const {
        auto it = edgeWeights.find({source, destination});
        if(it != edgeWeights.end()){
            return it->second;
        }
        return 0;
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
    // Function to get the neighbors of a vertex
    std::vector<int> getNeighborsKL(const int& vertexID) const{
        const auto it = adjacencyList.find(vertexID);
        if(it != adjacencyList.end()){
            const auto value = it->second;
            return value;
        }
        return {};

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

        std::cout << "AdjacencyLists:" << std::endl;
        for (const auto& entry : adjacencyList) {
            int key = entry.first;
            const std::vector<int>& neighbors = entry.second;

            std::cout << "Vertex " << key << " neighbors: ";
            for (int neighbor : neighbors) {
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }
    }

    // Function to clear the graph and reset it to its initial state
    void clear() {
        // Clear the vertices, adjacency list, and edge weights
        vertices.clear();
        adjacencyList.clear();
        edgeWeights.clear();

        // Clear mappings and other data
        coarserToFinerMappings.clear();
        edgesMappings.clear();
        verticesWeightsMapping.clear();
    }

    std::vector<int> inflateVertex(const int& vertexID, int& offset) {
        if(vertexID < 0) return std::vector<int>({vertexID});
        std::vector<int> inflation;
        auto vertex_it = vertices.find(vertexID);
        if(vertex_it == vertices.end()) exit(1);
        int weight = static_cast<int>(vertex_it->second);
        inflation.push_back(vertex_it->first);
        for(int i = 0; i < weight-1; i++){
            int new_id = -(offset+1);
            offset++;
            inflation.push_back(new_id);
            addVertex(new_id);
        }
        for(auto source: inflation){
            for(auto destination: inflation){
                if(destination != source){
                    edgeWeights[{source, destination}] = INF_EDGE_WEIGHT;
                    adjacencyList[source].push_back(destination);
                }
            }
        }
        return inflation;
    }
    void deflate(){
        for(auto it = vertices.begin(); it != vertices.end();){
            if(it->first < 0){
                it = vertices.erase(it);
            } else {
                ++it;
            }
        }
        for(auto it = adjacencyList.begin(); it != adjacencyList.end();){
            if(it->first < 0){
                it = adjacencyList.erase(it);
            } else {
                if(it->second.size()){
                    it->second.erase(std::remove_if(it->second.begin(), it->second.end(), [](int n) { return n < 0; }), it->second.end());
                }
                ++it;
            }
        }
        for(auto it = edgeWeights.begin(); it != edgeWeights.end();){
            if(it->second < 0){
                it = edgeWeights.erase(it);
            } else {
                ++it;
            }
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