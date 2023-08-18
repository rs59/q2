#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>

int coarsest_graph_size = 5; //TO BE DEFINED

// Function to compute matching between vertices for coarsening
std::unordered_map<int, int> ComputeMatching(Graph& graph_cm) {
    std::unordered_map<int, int> matching;
    std::unordered_set<int> unmatchedVertices;
    std::unordered_map<int, double> vertices = graph_cm.getVertices();

    // Populate the set of unmatched vertices
    for (const auto& vertex : vertices) {
        unmatchedVertices.insert(vertex.first);
    }

    for (const auto& vertex : vertices) {
        int vertexId = vertex.first;
        if (matching.find(vertexId) != matching.end()) {
            continue; // Skip already matched vertices
        }

        const std::vector<int>& neighbors = graph_cm.getNeighbors(vertexId);
        int matched_neighbor = -1;

        // Find the first unmatched neighbor
        for (int neighbor : neighbors) {
            if (unmatchedVertices.find(neighbor) != unmatchedVertices.end()) {
                matched_neighbor = neighbor;
                break;
            }
        }

        if (matched_neighbor != -1) {
            matching[vertexId] = matched_neighbor;
            matching[matched_neighbor] = vertexId;
            unmatchedVertices.erase(vertexId);
            unmatchedVertices.erase(matched_neighbor);
        }
    }

    return matching;
}


Graph CollapseVertices(Graph& graph_cv, std::unordered_map<int, int>& matching) {
    Graph coarsed_graph;

    //save in the new coarsed_graph the previous coarsening mappings
    for(int i = 0; i < graph_cv.getCoarsingLevel(); i++){
        coarsed_graph.pushBackMapping(graph_cv.getMapping(i));
        coarsed_graph.pushBackMappingEdges(graph_cv.getMappingEdges(i));
        coarsed_graph.pushBackVerticesWeights(graph_cv.getMappingVerticesWeights(i));
    }

    std::unordered_map<int,int> vertices_map;  //store the mapping between coarsed vertices and original ones {original vertex, new vertex}
    std::unordered_map<int, double> vertices = graph_cv.getVertices();

    int i = 0; //the new vertices numeration for coarsed graph
    for (const auto& vertex : vertices) {
        int vertexId = vertex.first;
        if (vertices_map.find(vertexId) != vertices_map.end()) {
            continue; // Skip already inserted vertices
        }
        auto matched_vertex_it = matching.find(vertexId);

        if (matched_vertex_it != matching.end()) {
            // Key (vertexID) was found in the map
            int matched_neighbor = matched_vertex_it->second; // Get the matched neighbor
            double new_weight = graph_cv.getVertexWeight(vertexId) + graph_cv.getVertexWeight(matched_neighbor);
            coarsed_graph.addVertex(i, new_weight);  //collapse the vertices
            vertices_map[vertexId] = i;
            vertices_map[matched_neighbor] = i;
            i = i + 1;
        } else {
            // Key (vertexID) was not found in the map
            coarsed_graph.addVertex(i, graph_cv.getVertexWeight(vertexId));
            vertices_map[vertexId] = i;
            i = i + 1;
        }
    }

    coarsed_graph.pushBackMapping(vertices_map);
    coarsed_graph.pushBackMappingEdges(graph_cv.getEdgeWeights());
    coarsed_graph.pushBackVerticesWeights(graph_cv.getVertices());

    return coarsed_graph;
}

// Function to update the edge weights after collapsing vertices
void UpdateEdgeWeights(Graph& graph_ue, const std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights) {
    std::unordered_map<int,int> vertices_map = graph_ue.getMapping(graph_ue.getCoarsingLevel() - 1);

    std::unordered_map<std::pair<int, int>, double, HashPair> edge_weights_updated;

    for (const auto& entry : edge_weights) {
        int nodeA = entry.first.first;
        int nodeB = entry.first.second;
        double value = entry.second;

        if (vertices_map.find(nodeA) != vertices_map.end() && vertices_map.find(nodeB) != vertices_map.end()) {
            int coarsedSource = vertices_map.at(nodeA);
            int coarsedDestination = vertices_map.at(nodeB);

            // Skip processing if coarsedSource and coarsedDestination are the same
            if (coarsedSource == coarsedDestination) {
                continue;
            }

            // Check if an edge already exists between the coarsed vertices
            if (graph_ue.containsEdge(coarsedSource, coarsedDestination)) {
                double existingWeight = graph_ue.getEdgeWeight(coarsedSource, coarsedDestination);
                // Sum the weights if an edge already exists
                graph_ue.addEdge(coarsedSource, coarsedDestination, existingWeight + (value/2)); //value is /2 because we have double edges
            } else {
                graph_ue.addEdge(coarsedSource, coarsedDestination, (value/2));
            }
        }
    }
}

Graph Coarsening(Graph& graph){
    Graph coarsened_graph = graph;

    while(coarsened_graph.size() > coarsest_graph_size){
        Graph temp_graph = coarsened_graph;

        std::unordered_map<int, int> matching = ComputeMatching(temp_graph);

        coarsened_graph = CollapseVertices(temp_graph, matching);   //graph with collapsed vertices and no edges

        std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights = temp_graph.getEdgeWeights();

        UpdateEdgeWeights(coarsened_graph, edgeWeights);
    }

    return coarsened_graph;

}