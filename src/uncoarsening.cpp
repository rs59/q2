#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>

Graph Uncoarsening(Graph& graph){
    int num_levels = graph.getCoarsingLevel();
    Graph uncoarsened_graph;

    for (int level = num_levels - 1; level >= 0; level--) {
        std::unordered_map<int, int> coarser_to_finer_mapping = graph.getMapping(level);
        std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping = graph.getMappingEdges(level);
        std::unordered_map<int, double> verticesWeights = graph.getMappingVerticesWeights(level);

        // Create a new graph to store the uncoarsened version
        Graph uncoarsened_temp;

        for (const auto &mapping: coarser_to_finer_mapping) {
            int original = mapping.first;

            auto vertexIt = verticesWeights.find(original);
            if(vertexIt != verticesWeights.end()){
                double weight = vertexIt->second;
                uncoarsened_temp.addVertex(original, weight);
            }
        }

        for (const auto& entry : edgesMapping) {
            int nodeA = entry.first.first;
            int nodeB = entry.first.second;
            double value = entry.second;

            uncoarsened_temp.addEdge(nodeA, nodeB, value);
        }

        // Update the uncoarsened graph for this level
        uncoarsened_graph = uncoarsened_temp;
    }

    return uncoarsened_graph;
}
