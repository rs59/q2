// Enable or disable debug output
#define DEBUG 1

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif


#include <iostream>
#include "KLCore.cpp"
#include "Reader.cpp"
#include <chrono>
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

void printPartitions(const std::vector<std::vector<int>>& partitions, std::unordered_map<int, double> vertices){
    int weight = 0;
    for(int i = 0; i< partitions.size(); i++){
        std::cout << "Partition " << i+1 << ": " << partitions[i].size() << std::endl;
        for(const auto& node: partitions[i]){
            std::cout << node << " ";
            weight += vertices[node];
        }
        std::cout << " weight: " << weight;
        std::cout << std::endl;

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

std::vector<std::vector<int>> makeNodePartion(const Graph& G, const int& numPartitions){
    auto nodes = G.getVertices();
    std::vector<std::vector<int>> partitions(numPartitions);
    int i = 0;
    for(auto& node : nodes){
        partitions[i].push_back(node.first);
        i = (i+1)%numPartitions;
    }
    /*i = 0;
    for(unsigned int j = 1; j< partitions.size(); j++){
        if(partitions[j].size() < partitions[0].size()){
            partitions[j].push_back(-(++i));
        }
    }*/
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
    return partitions;
}

std::vector<std::vector<int>> makePartion(const Graph& G, const int& numPartitions){
    std::vector<std::vector<int>> fPartitions(numPartitions);
    std::vector<Partition> partitions(numPartitions);
    double totVertexWeight = 0;
    auto vertices = G.getVertices();
    std::cout << "Number of Vertices :" << vertices.size()<< std::endl;
    /*for(auto vertex: vertices)
        std::cout << vertex.first << "," <<vertex.second << "  ";
    std::cout << std::endl;*/
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
    /*
    i = 0; 
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
    std::cout << "Initial CutSize: " << calculateCutSize(G, cleanPartition(partitions)) << std::endl;
    double max_weight = 0;
    for(auto& partition: partitions)
        max_weight = (partition.weight > max_weight) ? partition.weight : max_weight;
    int inserted = 0;
    for(auto& partition: partitions)
        fillPartition(partition, inserted, max_weight);
    /*
    for(int j = 0; j < numPartitions; j++) {
        for(const auto& node: partitions[j].part){
            std::vector<int> inflatedNodes = G.inflateVertex(node, inserted);
            fPartitions[j].insert(fPartitions[j].end(), inflatedNodes.begin(), inflatedNodes.end());
        }
    }
    return fPartitions;*/
    return cleanPartition(partitions);
}

void deflate(Graph& G,std::vector<std::vector<int>>& partitions){
    G.deflate();
    for(auto& partition: partitions)
        partition.erase(std::remove_if(partition.begin(), partition.end(), [](int n) { return n < 0; }), partition.end());
}

std::vector<std::vector<int>> multipartitionKL(Graph& G, const int& numPartitions) {
    auto partitions = makeNodePartion(G, numPartitions);
    //printPartitions(partitions);
    //if(1) return partitions;
    char optimized = 1;
    int rounds = 0;
    do {
        rounds++;
        //std::cout << "Round " << rounds <<std::endl;
        optimized = 1;
        for(int i = 0; i < numPartitions; i++){
            for(int j = i+1; j < numPartitions; j++){
                char inequal_size = 0;
                if(partitions[i].size() < partitions[j].size()){
                    partitions[i].push_back(-1);
                    inequal_size = 1;
                } else if (partitions[i].size() > partitions[j].size()){
                    partitions[j].push_back(-1);
                    inequal_size = 1;
                }
                int iterations = KL_Partitioning(G, partitions[i], partitions[j]);
                if(inequal_size){
                    partitions[i].erase(std::remove(partitions[i].begin(), partitions[i].end(), -1), partitions[i].end());
                    partitions[j].erase(std::remove(partitions[j].begin(), partitions[j].end(), -1), partitions[j].end());
                }
                if(iterations > 1 && optimized) optimized = 0;
            }
        }
    } while (!optimized);
    std::cout << "Number of Rounds: " << rounds << std::endl;
    //deflate(G, partitions);
    return partitions;
}

int main() {
    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    //resources\metismodels\x100000y220000m20q20.metis
    const string filename = "/content/q2/resources/metismodels/x10000y20000m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    //std::cout << "Registred Number of Vertices: " << graph.numVertices()<< std::endl;
    
    auto optPartitions = multipartitionKL(graph, NPART);
    printPartitions(optPartitions, graph.getVertices());
    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Execution time: " << seconds << " seconds" << std::endl;
    std::cout << "Final CutSize: " << calculateCutSize(graph, optPartitions);
    return 0;
}