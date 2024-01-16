#include "graph.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#define N 3

// Structure to represent a pair of nodes and the weight of the edge between them
struct SwapNodes {
    std::pair<int, int> nodes;
    double weight;
};

// Function to calculate the D values for a partition
std::map<int, double> DCalcPartition(const Graph &G, const std::vector<int> &nodesA, const std::vector<int> &nodesB) {
    std::map<int, double> dValues;
    // Iterate through nodes in nodesA
    for (auto &nodeA : nodesA) {
        double tDValue = 0;
        auto nNodes = G.getNeighborsKL(nodeA);
        // Iterate through neighbors of the current node
        for (auto node : nNodes) {
            auto tEdgeWeight = G.getEdgeWeightKL(nodeA, node);
            auto nodeAIt = std::find(nodesA.begin(), nodesA.end(), node);
            auto nodeBIt = std::find(nodesB.begin(), nodesB.end(), node);
            // Update the D value based on the partition of the neighboring node
            if (nodeAIt != nodesA.end()) {
                tDValue -= tEdgeWeight;
            } else if (nodeBIt != nodesB.end()) {
                tDValue += tEdgeWeight;
            }
        }
        dValues.insert(std::pair<int, double>(nodeA, tDValue));
    }
    return dValues;
}

// Function to calculate D values for both partitions
void calcD(const Graph &G, const std::vector<int> &nodesA, std::map<int, double> &dValuesA, const std::vector<int> &nodesB, std::map<int, double> &dValuesB) {
    dValuesA = DCalcPartition(G, nodesA, nodesB);
    dValuesB = DCalcPartition(G, nodesB, nodesA);
}

// Function to update D values after a node swap
void updateD(const Graph &G, std::map<int, double> &dValuesA, std::map<int, double> &dValuesB, std::pair<int, int> swapNodes) {
    dValuesA.erase(swapNodes.first);
    dValuesB.erase(swapNodes.second);
    for (auto it = dValuesA.begin(); it != dValuesA.end(); ++it) {
        it->second += 2 * G.getEdgeWeightKL(swapNodes.first, it->first);
        it->second -= 2 * G.getEdgeWeightKL(swapNodes.second, it->first);
    }
    for (auto it = dValuesB.begin(); it != dValuesB.end(); ++it) {
        it->second += 2 * G.getEdgeWeightKL(swapNodes.second, it->first);
        it->second -= 2 * G.getEdgeWeightKL(swapNodes.first, it->first);
    }
}

// Function to calculate gains for all possible node swaps
std::multimap<double, std::pair<int, int>> gainCalc(Graph &G, std::map<int, double> dValuesA, std::map<int, double> dValuesB) {
    std::multimap<double, std::pair<int, int>> gains;
    for (auto dValueA : dValuesA) {
        for (auto dValueB : dValuesB) {
            double gain = dValueA.second + dValueB.second - 2 * G.getEdgeWeightKL(dValueA.first, dValueB.first);
            gains.insert(std::pair<double, std::pair<int, int>>(gain, std::pair<int, int>(dValueA.first, dValueB.first)));
        }
    }
    return gains;
}

// Function to calculate the index with maximum gain in a list of swap nodes
int calcMaxGain(std::vector<SwapNodes> swapNodes, double &max_gain) {
    double gain;
    int k;
    for (int i = 0; i < swapNodes.size(); i++) {
        if (i == 0) {
            k = i;
            gain = swapNodes[i].weight;
            max_gain = swapNodes[i].weight;
        } else {
            gain += swapNodes[i].weight;
            if (gain > max_gain) {
                max_gain = gain;
                k = i;
            }
        }
    }
    return k;
}

// Function to print gains
void printGain(std::multimap<double, std::pair<int, int>> gains) {
    // std::cout << "Gains: " << std::endl;
    // for (auto gain : gains) {
    //   if (gain.first != 0 )
        // std::cout << "Edge (" << gain.second.first << "," << gain.second.second << "): " << gain.first << std::endl;
    // }
}

// Function to update the partition after a swap
void updatePartition(std::vector<int> &nodesA, std::vector<int> &nodesB, std::vector<SwapNodes> swapNodes, int k) {
    for (int i = 0; i < k + 1; i++) {
        auto itA = std::find(nodesA.begin(), nodesA.end(), swapNodes[i].nodes.first);
        auto itB = std::find(nodesB.begin(), nodesB.end(), swapNodes[i].nodes.second);
        if (itA != nodesA.end() && itB != nodesB.end()) {
            std::iter_swap(itA, itB);
        }
    }
}

// Function to print the partition
void partitionPrint(std::vector<int> nodesA, std::vector<int> nodesB) {
    std::cout << "Partition A: " << std::endl;
    for (auto i : nodesA) {
        std::cout << i << ", ";
    }
    std::cout << std::endl
              << "Partition B: " << std::endl;
    for (auto i : nodesB) {
        std::cout << i << ", ";
    }
    std::cout << std::endl;
}

// Main KL Partitioning function
int KL_Partitioning(Graph &G, std::vector<int> &nodesA, std::vector<int> &nodesB) {
    std::map<int, double> dValuesA, dValuesB;
    double max_gain;
    int i = 0;
    // Repeat until convergence or maximum iterations
    do {
        i++;
        std::vector<SwapNodes> swapNodes;
        calcD(G, nodesA, dValuesA, nodesB, dValuesB);
        // Repeat until no further improvement in gain
        while (!dValuesA.empty() && !dValuesB.empty()) {
            auto gains = gainCalc(G, dValuesA, dValuesB);
            auto mostGain_it = gains.rbegin();
            if (mostGain_it == gains.rend()) {
                exit(-1);
            }
            swapNodes.push_back({mostGain_it->second, mostGain_it->first});
            updateD(G, dValuesA, dValuesB, mostGain_it->second);
            // std::cout << "\tc  " << std::endl;
            printGain(gains);
        }
        // Find the index with the maximum cumulative gain
        int k = calcMaxGain(swapNodes, max_gain);
        std::cout << "\tmg " << max_gain << std::endl;
        // Update partition if there is a positive gain
        if (max_gain > 0)
            updatePartition(nodesA, nodesB, swapNodes, k);
    } while (max_gain > 0 && i < 100); // Convergence condition
    return i; // Return the number of iterations
}


/*
using namespace std;

void printPart(vector<int> nodesA, vector<int> nodesB) {
    cout << "Partition A :" << endl;
    for(auto node: nodesA){
        cout << node << " ";
    }
    cout << endl;
    cout << "Partition B :" << endl;
    for(auto node: nodesB){
        cout << node << " ";
    }
    cout << endl;
}


int main()
{
    Graph G;
    G.addVertex(1, 1);
    G.addVertex(2, 1);
    G.addVertex(3, 1);
    G.addVertex(4, 1);
    G.addVertex(5, 1);
    G.addVertex(6, 1);

    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 3);
    G.addEdge(1, 5, 2);
    G.addEdge(1, 6, 4);

    G.addEdge(2, 3, 1);
    G.addEdge(2, 4, 4);
    G.addEdge(2, 5, 2);
    G.addEdge(2, 6, 1);

    G.addEdge(3, 4, 3);
    G.addEdge(3, 5, 2);
    G.addEdge(3, 6, 1);

    G.addEdge(4, 5, 4);
    G.addEdge(4, 6, 3);

    G.addEdge(5, 6, 2);
    vector<int> nodesA = {1, 2, 3}, nodesB = {4, 5, 6};
    printPart(nodesA,nodesB);
    auto iter = KL_Partitioning(G, nodesA, nodesB);
    cout << endl << "iterations: " << iter << endl;
    printPart(nodesA,nodesB);
    iter = KL_Partitioning(G, nodesA, nodesB);
    cout << endl << "iterations: " << iter << endl;
    printPart(nodesA,nodesB);
    return 0;
}**/