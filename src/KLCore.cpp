#ifdef DEBUG
#include <iostream>
#endif
#include "graph.h"
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
    // add or subtract based on the partition the value of the edge
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
            // Formula for the gain
            double gain = dValueA.second + dValueB.second - 2 * G.getEdgeWeightKL(dValueA.first, dValueB.first);
            // std::cout << "gain " << gain << " from " << dValueA.first << " and " << dValueB.first << std::endl;
            // save couple of nodes and their gain in ordered map
            gains.insert(std::pair<double, std::pair<int, int>>(gain, std::pair<int, int>(dValueA.first, dValueB.first)));
        }
    }
    return gains;
}

// Function to calculate the index with maximum gain in a list of swap nodes
int calcMaxGain(std::vector<SwapNodes> swapNodes, double &max_gain) {
    double gain;
    int k;
    // sum all single max_gain and save the partial sum with the highest sum
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

#ifdef DEBUG
// Function to print gains
void printGain(std::multimap<double, std::pair<int, int>> gains) {
    std::cout << "Gains: " << std::endl;
    for (auto gain : gains) {
    if (gain.first != 0 )
        std::cout << "Edge (" << gain.second.first << "," << gain.second.second << "): " << gain.first << std::endl;
    }
}
#endif

// Function to update the partition after a swap
void updatePartition(std::vector<int> &nodesA, std::vector<int> &nodesB, std::vector<SwapNodes> swapNodes, int k) {
    for (int i = 0; i < k + 1; i++) {
        // for all memebers which gave the max amount of swap find the iterator
        auto itA = std::find(nodesA.begin(), nodesA.end(), swapNodes[i].nodes.first);
        auto itB = std::find(nodesB.begin(), nodesB.end(), swapNodes[i].nodes.second);
        if (itA != nodesA.end() && itB != nodesB.end()) {
            // swap the nodes in the partitions
            std::iter_swap(itA, itB);
        }
    }
}

#ifdef DEBUG
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
#endif

// Main KL Partitioning function
int KL_Partitioning(Graph &G, std::vector<int> &nodesA, std::vector<int> &nodesB) {
    std::map<int, double> dValuesA, dValuesB;
    double max_gain;
    int i = 0;
    // Repeat until convergence or maximum iterations
    do {
        i++;
        std::vector<SwapNodes> swapNodes;
        DEBUG_STDOUT("calcD ");
        calcD(G, nodesA, dValuesA, nodesB, dValuesB);
        DEBUG_STDOUT("finished calcD ");
        int prevPrevGain = -3000000;
        int prevGain = -2000000;
        int currGain = -1000000;
        // Repeat until no further improvement in gain
        while (!dValuesA.empty() && !dValuesB.empty()) {
            // std::cout << "gainCalc " << std::endl;
            auto gains = gainCalc(G, dValuesA, dValuesB);
            // std::cout << "finished gainCalc " << std::endl;
            auto mostGain_it = gains.rbegin();
            if (mostGain_it == gains.rend()) {
                exit(-1);
            }
            swapNodes.push_back({mostGain_it->second, mostGain_it->first});
            updateD(G, dValuesA, dValuesB, mostGain_it->second);
            // std::cout << "\tc  " << std::endl;
            // prevPrevGain = prevGain;
            // prevGain = currGain;
            // currGain = (*gains.rend()).first;
            // std::cout << "prevPrevgain " << prevPrevGain << "prevgain " << prevGain << "Currgain " << currGain << std::endl;
            // if(prevGain==prevPrevGain && prevGain==currGain) {
            //   std::cout << "repeated gains, continuing  " << std::endl;
            //   break;
            // }
        }
        // Find the index with the maximum cumulative gain
        int k = calcMaxGain(swapNodes, max_gain);
        DEBUG_STDOUT("\tMax Gain: " +  std::to_string(max_gain));
        // std::cout << "\tmg " << max_gain << std::endl;
        // Update partition if there is a positive gain
        if (max_gain > 0)
            updatePartition(nodesA, nodesB, swapNodes, k);
    } while (max_gain > 0 && i < 100); // Convergence condition
    return i; // Return the number of iterations
}

