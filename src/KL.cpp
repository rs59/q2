#include "graph copy.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

struct Partition {
    std::vector<int> nodes;
    double weight;
};

struct SwapNodes {
    std::pair<int, int> nodes;
    double weight;
};



std::map<int, double> DCalcPartition(const Graph &G, const Partition& a, const Partition& b) {
    std::map<int, double> dValues;
    for (auto &nodeA : a.nodes) {
        double dValue = 0;
        auto nNodes = G.getNeighbors(nodeA);
        for (auto node : nNodes) {
            // The Edge Exist as i'm exploring only the neibourghs
            auto tEdgeWeight = G.getEdgeWeight(nodeA, node);
            auto nodeAIt = std::find(a.nodes.begin(), a.nodes.end(), node);
            auto nodeBIt = std::find(b.nodes.begin(), b.nodes.end(), node);
            if(nodeAIt != a.nodes.end()) {
                    dValue -= tEdgeWeight;
            }else if(nodeBIt != b.nodes.end()) {
                    dValue += tEdgeWeight;
            }
        }
        dValues.insert(std::pair<int, double>(nodeA, dValue));
    }
    return dValues;
}

void calcD(const Graph& G, const Partition& a, const Partition& b, std::map<int, double>& dValuesA, std::map<int, double>& dValuesB) {
    dValuesA = DCalcPartition(G, a, b);
    dValuesB = DCalcPartition(G, b, a);
}

void updateD(const Graph &G, std::map<int, double> &dValuesA, std::map<int, double> &dValuesB, std::pair<int, int> swapNodes) {
    dValuesA.erase(swapNodes.first);
    dValuesB.erase(swapNodes.second);
    for (auto it = dValuesA.begin(); it != dValuesA.end(); ++it) {
        it->second += 2 * G.getEdgeWeight(swapNodes.first, it->first);
        it->second -= 2 * G.getEdgeWeight(swapNodes.second, it->first);
    }
    for (auto it = dValuesB.begin(); it != dValuesB.end(); ++it) {
        it->second += 2 * G.getEdgeWeight(swapNodes.second, it->first);
        it->second -= 2 * G.getEdgeWeight(swapNodes.first, it->first);
    }
}

void dValPrint(std::map<int, double> dValuesA, std::map<int, double> dValuesB)
{
    std::cout << "D values A: " << std::endl;
    for (auto dVal : dValuesA)
    {
        std::cout << dVal.first << ": " << dVal.second << std::endl;
    }
    std::cout << std::endl << "D values B: " << std::endl;
    for (auto dVal : dValuesB)
    {
        std::cout << dVal.first << ": " << dVal.second << std::endl;
    }
    std::cout << std::endl;
}

// }
std::multimap<double, std::pair<int, int>> gainCalc(Graph &G, std::map<int, double> dValuesA, std::map<int, double> dValuesB) {
    std::multimap<double, std::pair<int, int>> gains;
    for (auto dValueA : dValuesA) {
        for (auto dValueB : dValuesB) {
            double gain = dValueA.second + dValueB.second - 2 * G.getEdgeWeight(dValueA.first, dValueB.first);
            gains.insert(std::pair<double, std::pair<int, int>>(gain, std::pair<int, int>(dValueA.first, dValueB.first)));
        }
    }
    return gains;
}

int calcMaxGain(std::vector<SwapNodes> swapNodes, double &max_gain) {
    double gain = swapNodes[0].weight;;
    int k = 0;
    max_gain = swapNodes[0].weight;
    for (int i = 1; i < swapNodes.size(); i++) {
        gain += swapNodes[i].weight;
        if (gain > max_gain) {
            max_gain = gain;
            k = i;
        }
    }
    return k;
}

void printGain(std::multimap<double, std::pair<int, int>> gains){
    std::cout << "Gains: " << std::endl; 
    for(auto gain: gains){
        std::cout << "Edge ("<< gain.second.first << "," << gain.second.second << "): " << gain.first << std::endl;
    }
}

void updatePartition(Partition& a, Partition& b, std::vector<SwapNodes>& swapNodes, int k) {
    for (int i = 0; i <= k; i++) {
        // Find the indices of the values in the vectors
        auto itA = std::find(a.nodes.begin(), a.nodes.end(), swapNodes[i].nodes.first);
        auto itB = std::find(b.nodes.begin(), b.nodes.end(), swapNodes[i].nodes.second);

        // Check if the values were found
        if (itA != a.nodes.end() && itB != b.nodes.end()) {
            // Swap the values using std::swap
            std::iter_swap(itA, itB);
        } else {
            exit(-2);
        }
    }
}

void partitionPrint(const Partition& a, const Partition& b){
    std::cout << "Partiotion A: " <<std::endl;
    for(auto i: a.nodes){
        std::cout << i << ", ";
    }
    std::cout << std::endl << "Partiotion B: " <<std::endl;
    for(auto i: b.nodes){
        std::cout << i << ", ";
    }
    std::cout << std::endl;
}


int KL_Partitioning(Graph &G, Partition& a, Partition& b) {
    std::map<int, double> dValuesA, dValuesB;
    double max_gain;
    int i = 0;
    do {
        i++;
        std::vector<SwapNodes> swapNodes;
        calcD(G, a, b, dValuesA, dValuesB);
        while (!dValuesA.empty() && !dValuesB.empty()) {
            auto gains = gainCalc(G, dValuesA, dValuesB);
            auto mostGain_it = gains.rbegin();
            if (mostGain_it == gains.rend()) {
                exit(-1);
            }
            swapNodes.push_back({mostGain_it->second, mostGain_it->first});
            updateD(G, dValuesA, dValuesB, mostGain_it->second);
        }
        int k = calcMaxGain(swapNodes, max_gain);
        std::cout << max_gain << ", " << k << " ";
        if(max_gain > 0)
            updatePartition(a, b, swapNodes, k);
        swapNodes.clear();
    } while (max_gain > 0);
    std::cout << std::endl;
    return i;
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
