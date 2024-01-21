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
        expandedRange = other.expandedRange;
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
        expandedRange = other.expandedRange;
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

    Graph& operator=(Graph&& other) noexcept;

    // Function to add a vertex to the graph with its weight
    void addVertex(int vertexID, double weight = 1) {
        vertices[vertexID] = weight;
        adjacencyList[vertexID] = std::vector<int>();

    }

    void setOriginalVertices(int originalVertices) {
        // Update ORIGINAL_VERTICES_COUNT
        ORIGINAL_VERTICES_COUNT = originalVertices;
    }

    
    // Function to add a vertex to the graph with its weight
    void addExpandedVertex(int parentVertexID, double weight = 1) {
        int thisID = ORIGINAL_VERTICES_COUNT;
        vertices[thisID] = weight;
        adjacencyList[thisID] = std::vector<int>();
    }

    int getExpandedStart() const {
      return ORIGINAL_VERTICES_COUNT;
    }

    // Expand all nodes: pass their weights out to their children
    void expandNodes() {
      int startingPosition = getExpandedStart() + 1;
      // std::cout << "starting pos: " << startingPosition << std::endl;
      auto vertices_unmod = vertices;
      for(const auto& entry : vertices_unmod){ // loop through all nodes
        int parentNodeID = entry.first;
        double weight = entry.second;
        // std::cout << "ID: " << parentNodeID << "weight: " << weight << std::endl;
        if(weight >= 1) {
          expandedRange[parentNodeID]=startingPosition;
          for(int i=0; i < weight; i++) {  // if it is an expanded node
            vertices[startingPosition+i]=1; // Step 2: Assign all extra generated vertices to the same partition as their original vertex
            // std::cout << "new node id: " << startingPosition+i << std::endl;
            expandedRange[startingPosition+i]=parentNodeID;
            addEdge((i == 0) ? parentNodeID: startingPosition+i-1, startingPosition+i, static_cast<double>(40));
          }
          startingPosition += weight;
        } else {
          expandedRange[parentNodeID]=-1; // mark non-expanded nodes by -1
        }
      }
    }
    
    // Remove all nodes: pass their weights out to their children
    void collapseNodes() {
      
      for(const auto& entry : vertices){ // loop through all nodes
        int expandedNodeID = entry.first;
        if(expandedNodeID>=getExpandedStart()) { // then it is an expanded node and can be deleted
          vertices.erase(expandedNodeID);
          adjacencyList.erase(expandedNodeID);
        }
      }

    }

    // Function to get a vector of unique vertex IDs
    // std::vector<int> getUniqueVertices() const
    // {
    //     std::vector<int> uniqueVertices;

    //     for (const auto &entry : vertices)
    //     {
    //         uniqueVertices.push_back(entry.first);
    //     }

    //     return uniqueVertices;
    // }

    // Function to split the graph into two based on given vectors of vertices
    std::pair<Graph, Graph> splitGraph(const std::vector<int> &vertices1, const std::vector<int> &vertices2) 
    {
        // Create two new graphs
        Graph graph1, graph2;

        // Copy vertices and weights to the new graphs
        for (int vertexID : vertices1)
        {
            if (containsVertex(vertexID))
            {
                double weight = getVertexWeight(vertexID);
                graph1.addVertex(vertexID, weight);

                // Copy neighbors and edges
                for (int neighbor : getNeighbors(vertexID))
                {
                    if (std::find(vertices1.begin(), vertices1.end(), neighbor) != vertices1.end())
                    {
                        // If the neighbor is in vertices1, add the edge to graph1
                        double edgeWeight = getEdgeWeight(vertexID, neighbor);
                        graph1.addEdge(vertexID, neighbor, edgeWeight);
                    }
                }
            }
        }

        for (int vertexID : vertices2)
        {
            if (containsVertex(vertexID))
            {
                double weight = getVertexWeight(vertexID);
                graph2.addVertex(vertexID, weight);

                // Copy neighbors and edges
                for (int neighbor : getNeighbors(vertexID))
                {
                    if (std::find(vertices2.begin(), vertices2.end(), neighbor) != vertices2.end())
                    {
                        // If the neighbor is in vertices2, add the edge to graph2
                        double edgeWeight = getEdgeWeight(vertexID, neighbor);
                        graph2.addEdge(vertexID, neighbor, edgeWeight);
                    }
                }
            }
        }

        // Add connections for disconnected nodes within their own half
        addDisconnectedConnections(graph1);
        addDisconnectedConnections(graph2);

        return {graph1, graph2};
    }

// Used to make a temp graph gain the mapping informations for uncoarsening thte new level
void copyCoarseningData(Graph &other, int level)
{
    // Clear mappings and other data
    coarserToFinerMappings.clear();
    edgesMappings.clear();
    verticesWeightsMapping.clear();

    int i=0;

    for (std::unordered_map<int, int> coarser_to_finer_mapping : other.coarserToFinerMappings)
    {
        if(i < level){
            this->pushBackMapping(coarser_to_finer_mapping);
        }
        i++;
    }
    i = 0;
    for (std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping : other.edgesMappings)
    {
        if(i < level){
            this->pushBackMappingEdges(edgesMapping);
        }
        i++;
    }
    i = 0;
    for (std::unordered_map<int, double> verticesWeights : other.verticesWeightsMapping)
    {
        if(i < level){
            this->pushBackVerticesWeights(verticesWeights);
        }
        i++;
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

    void setVertexWeight(const int& vertexID,const double& weight){
        vertices[vertexID] = weight;
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
    std::unordered_map<int, double> getVertices() const {
        return vertices;
    }

    //Function to get the edges list and weights
    std::unordered_map<std::pair<int, int>, double, HashPair> getEdgeWeights(){
        return edgeWeights;
    }

    std::unordered_map<int, int> getExpandedRange() const {
      return expandedRange;
    }

    //get the number of levels
    int getCoarsingLevel(){
        return coarserToFinerMappings.size();
    }

    //Insert vertices mapping for a new level
    void pushBackMapping(std::unordered_map<int,int>mapping){
        coarserToFinerMappings.push_back(mapping);
    }

    std::unordered_map<int,int> getMapping(int level){
        return coarserToFinerMappings[level];
    }

    //Insert edges mapping for a new level
    void pushBackMappingEdges(std::unordered_map<std::pair<int, int>, double, HashPair> edges){
        edgesMappings.push_back(edges);
    }

    std::unordered_map<std::pair<int, int>, double, HashPair> getMappingEdges(int level){
        return edgesMappings[level];
    }

    //Insert vertices weights for a new level
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

        
        std::cout << "ExpandedRange:" << std::endl;
        for (const auto& entry : expandedRange) {
            const auto& key = entry.first;
            int expandedRangeVal = entry.second;

            std::cout << "Vertex: " << key << ", ExpandedRange: " << expandedRangeVal << std::endl;

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
        expandedRange.clear();
        
        // Reset ORIGINAL_VERTICES_COUNT
        ORIGINAL_VERTICES_COUNT = 0;
    }

    // Function to clear the graph and reset it to its initial state
    void clearButMappings() {
        // Clear the vertices, adjacency list, and edge weights
        vertices.clear();
        adjacencyList.clear();
        edgeWeights.clear();

        // Clear mappings and other data
        expandedRange.clear();

        // Reset ORIGINAL_VERTICES_COUNT
        ORIGINAL_VERTICES_COUNT = 0;
    }

    std::vector<int> inflateVertex(const int& vertexID, int& offset) {
        if(vertexID < 0) return std::vector<int>();
        std::vector<int> inflation;
        auto vertex_it = vertices.find(vertexID);
        if(vertex_it == vertices.end()) exit(1);
        int weight = static_cast<int>(vertex_it->second);
        //inflation.push_back(vertex_it->first);
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
    // Member variable to track the count of original vertices
    int ORIGINAL_VERTICES_COUNT;

    void addDisconnectedConnections(Graph &graph) const
    {
        std::vector<int> disconnectedNodes;

        // Find disconnected nodes in the graph
        for (const auto &entry : graph.getVertices())
        {
            int vertexID = entry.first;
            if (graph.getNeighbors(vertexID).empty())
            {
                // std::cout << "empty " << vertexID << std::endl;
                disconnectedNodes.push_back(vertexID);
            }
        }

        // Connect disconnected nodes within their own subgraph
        for (int disconnectedNode : disconnectedNodes)
        {
            // std::cout << " disconnected node search " << std::endl;
            for (const auto &entry : graph.getVertices())
            {
                int otherNode = entry.first;
                if (disconnectedNode != otherNode)
                {
                    graph.addEdge(disconnectedNode, otherNode, 1.0); // You can adjust the edge weight as needed
                    // std::cout << "added disconnected node " << disconnectedNode << " to " << otherNode << std::endl;
                    break;
                }
            }
        }
    }

    // Data structures to store the graph data
    std::unordered_map<int, double> vertices;  // VertexID -> Weight
    std::unordered_map<int, int> expandedRange;  // VertexID -> Starting vertex ID of expanded vertices (or 0, indicating none); weight (number of expanded vertices) comes from vertices unordered_map
    std::unordered_map<int, std::vector<int>> adjacencyList;  // VertexID -> List of neighbors
    std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights; // (Source, Destination) -> Edge Weight

    // Container to hold vertex and edges mappings for each level
    std::vector<std::unordered_map<int, int>> coarserToFinerMappings;
    std::vector<std::unordered_map<std::pair<int, int>, double, HashPair>> edgesMappings;
    std::vector<std::unordered_map<int, double>> verticesWeightsMapping;
};

Graph& Graph::operator=(Graph&& other) noexcept {
    if (this != &other) {
        // Move assign vertices, adjacencyList, and edgeWeights
        vertices = std::move(other.vertices);
        adjacencyList = std::move(other.adjacencyList);
        edgeWeights = std::move(other.edgeWeights);

        // Move assign other data structures
        coarserToFinerMappings = std::move(other.coarserToFinerMappings);
        edgesMappings = std::move(other.edgesMappings);
        verticesWeightsMapping = std::move(other.verticesWeightsMapping);
        expandedRange = std::move(other.expandedRange);

        // Move assign scalar members
        ORIGINAL_VERTICES_COUNT = other.ORIGINAL_VERTICES_COUNT;
    }
    return *this;
}

#endif // GRAPH_H