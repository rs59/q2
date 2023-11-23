#include <iostream>
#include "KLCore.cpp"
#include "reader.cpp"
#define NPART 7

struct Partition {
    std::vector<int> part;
    double weight;
};

double calculateCutSize(const Graph& G, const std::vector<std::vector<int>>& partitions){
    double cutSize = 0;
    for(auto& partition : partitions)
        for(auto& node : partition){
            auto nNodes = G.getNeighborsKL(node);
            for(auto& nNode:nNodes)
                if(std::find(partition.begin(), partition.end(), nNode) == partition.end())
                    cutSize += G.getEdgeWeightKL(node, nNode);
        }
    return cutSize/2;
}

void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight){
    for(int i = 0; i < static_cast<int>(max_weight-partition.weight); i++){
        alreadyIns++;
        partition.part.push_back(-alreadyIns);
    }
}

std::vector<std::vector<int>> cleanPartition(const std::vector<Partition>& partitions){
    std::vector<std::vector<int>> cleanPartitions;
    for(int i = 0; i < partitions.size(); i++){
        cleanPartitions.push_back(partitions[i].part);
        /*std::cout << "Partition " << i << ":" <<std::endl; 
        for(auto n: cleanPartitions[i])
            std::cout << n << " ";
        std::cout<<std::endl;*/
    }
    return cleanPartitions;
}

std::vector<std::vector<int>> makePartion(Graph & G, const int& numPartitions){
    std::vector<std::vector<int>> fPartitions(numPartitions);
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
        partitions[i].part.push_back(it ->second);
        partitions[i].weight += it->first;
        i = (i+1)%numPartitions;
    }
    /*i = 0; 
    for(auto p: partitions){
        i++;
        std::cout << "Partition " << i << " : weight " << p.weight<< std::endl;
        for(auto n: p.part)
            std::cout << n << " ";
        std::cout << std::endl;
    }
    auto tpart = cleanPartition(partitions);
    for(auto p: tpart){
        i++;
        std::cout << "Partition " << i << " :" << std::endl;
        for(auto n: p)
            std::cout << n << " ";
        std::cout << std::endl;
    }*/
    std::cout << "Initial CutSize: " << calculateCutSize(G, cleanPartition(partitions));
    double max_weight = 0;
    for(auto& partition: partitions)
        max_weight = (partition.weight > max_weight) ? partition.weight : max_weight;
    int inserted = 0;
    for(auto& partition: partitions)
        fillPartition(partition, inserted, max_weight);

    for(int j = 0; j < numPartitions; j++) {
        for(const auto& node: partitions[j].part){
            std::vector<int> inflatedNodes = G.inflateVertex(node, inserted);
            fPartitions[j].insert(fPartitions[j].end(), inflatedNodes.begin(), inflatedNodes.end());
        }
    }
    return fPartitions;
}

void deflate(Graph& G,std::vector<std::vector<int>>& partitions){
    G.deflate();
    for(auto& partition: partitions)
        partition.erase(std::remove_if(partition.begin(), partition.end(), [](int n) { return n < 0; }), partition.end());
}

std::vector<std::vector<int>> multipartitionKL(Graph& G, const int& numPartitions) {
    auto partitions = makePartion(G, numPartitions);
    char optimized = 1;
    int rounds = 0;
    do {
        rounds++;
        std::cout << std::endl << "Round " << rounds <<std::endl;
        optimized = 1;
        for(int i = 0; i < numPartitions; i++){
            for(int j = i+1; j < numPartitions; j++){
                int iterations = KL_Partitioning(G, partitions[i], partitions[j]);
                    std::cout << iterations << " ";
                if(iterations >= 100 && optimized) optimized = 0;
            }
        }
    } while (!optimized);
    deflate(G, partitions);
    return partitions;
}

int main() {
    //resources\metismodels\x100000y220000m20q20.metis
    const string filename = "resources/metismodels/x100y200m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    auto optPartitions = multipartitionKL(graph, NPART);
    std::cout << "Final CutSize: " << calculateCutSize(graph, optPartitions);
    return 0;
}
