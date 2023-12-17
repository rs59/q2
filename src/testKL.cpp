#include <map>
#include <iostream>
#include <iomanip>
#include "graph.h"
#include "reader.cpp"
#define NPART 7

struct Partition {
    std::vector<std::pair<double, int>> part;
    double weight;
};

void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight){
    for(int i = 0; i < static_cast<int>(max_weight-partition.weight); i++){
        alreadyIns++;
        partition.part.push_back(std::pair<double, int>(0, -alreadyIns));
    }
}

std::vector<vector<int>> makePartion(Graph & G, int numPartitions){
    std::vector<vector<int>> fPartitions(numPartitions);
    std::vector<Partition> partitions(numPartitions);
    double totVertexWeight = 0;
    auto vertices = G.getVertices();
    std::multimap<double, int> nodes;
    for(auto& node: vertices){
        totVertexWeight += node.second;
        nodes.insert(std::pair<double, int>(node.second, node.first));
    }
    vertices.clear();
    double nodeSize = totVertexWeight/static_cast<double>(numPartitions);
    double excessSize = static_cast<int>(totVertexWeight)%numPartitions;
    size_t i = 0;
    for(auto it = nodes.rbegin(); it != nodes.rend(); ++it){
        if(partitions[i].weight >= nodeSize-0.5){
            size_t j = (i+1)%numPartitions;
            while(partitions[j].weight >= nodeSize-0.5 && i != j) {
                j = (j+1)%numPartitions;
            }
            i = j;
        }
        partitions[i].part.push_back(std::pair(it->first, it ->second));
        partitions[i].weight += it->first;
        i = (i+1)%numPartitions;
    }
    double max_weight = 0;
    for(auto& partition: partitions)
        max_weight = (partition.weight > max_weight) ? partition.weight : max_weight;
    int inserted = 0;
    for(auto& partition: partitions)
        fillPartition(partition, inserted, max_weight);

    for(int j = 0; j < numPartitions; j++) {
        for(const auto& node: partitions[j].part){
            std::vector<int> inflatedNodes = G.inflateVertex(node.second, inserted);
            fPartitions[j].insert(fPartitions[j].end(), inflatedNodes.begin(), inflatedNodes.end());
        }
    }
    return fPartitions;
}

using namespace std;

int main(){
    //resources\metismodels\x100000y220000m20q20.metis
    const string filename = "resources/metismodels/x100y200m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    auto partitions = makePartion(graph, NPART);
    size_t i = 0;
    double tot_weight = 0;
    for(auto partition: partitions){
        i++;
        cout << "Partition " << i << " weight = " << partition.size() << endl;
        tot_weight += partition.size();
    }
    cout << "Total Weight: " << tot_weight << endl;
    // for(auto element: partitions[12].part){
    //     cout << "Node " << element.second << ": " << element.first << endl;
    // }
    //graph.print();
    /*int offset = 0;
    auto infl = graph.inflateVertex(98,offset);
    for(auto vert: infl){
        cout << vert << ", ";
    }
    cout << endl;
    auto nnodes = graph.getNeighbors(98);
    for(auto vert: nnodes){
        cout << vert << ", ";
    }
    nnodes.push_back(98);
    cout << endl << "\t";
    for(int destination: nnodes)
        cout << destination << "\t";
    cout << endl;
    for(int source: nnodes){
        cout << source << "\t";
        for(int destination: nnodes)
            cout << graph.getEdgeWeightKL(source, destination) << "\t";
        cout << endl;
    }

    cout << endl;*/

}