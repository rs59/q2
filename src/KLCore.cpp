#include "graph.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

multimap<int, int> DCalcPartition(Graph& G, vector<int> nodesA, vector<int> nodesB){
    multimap<int, int> dValuesA;
    for(auto& nodeA: nodesA){
        int tDValue = 0;
        auto nNodes = G.getNeighbors(nodeA);
        for(auto node:nNodes){
            //The Edge Exist as i'm exploring only the neibourghs
            auto tEdgeWeight = G.getEdgeWeight(nodeA, node);
            auto nodeIt = std::find(nodesB.begin(), nodesB.end(), node);
            if(nodeIt != nodesB.end()){
                tDValue += tEdgeWeight; 
            }
            else{
                tDValue -= tEdgeWeight;
            }
            
            //Can't calculate the B partition values cause of missing internal edges
        }
        dValuesA.insert(pair<int, int>(tDValue, nodeA));
    }
    return dValuesA;
}

void DCalc(Graph& G, vector<int> nodesA, vector<int> nodesB){
    auto dValuesA = DCalcPartition(G, nodesA, nodesB);
    auto dValuesB = DCalcPartition(G, nodesB, nodesA);
}


int main () {
    return 0;

}