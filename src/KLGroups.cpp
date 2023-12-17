#include <iostream>
#include "KL.cpp"

void printPartitions(const std::vector<Partition>& partitions){
    for(int i = 0; i< partitions.size(); i++){
        std::cout << "Partition " << i+1 << ": Nodes " << partitions[i].nodes.size() << ", Weight "<< partitions[i].weight << std::endl;
        for(const auto& node: partitions[i].nodes){
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}





/* 
 * Explore each node of the graph @param G
 * get for each of them all the neighbours: those with a connecting edge
 * The cutsize is the sum of all external edges referring to the partitions
 */
double calculateCutSize(const Graph& G, const std::vector<Partition>& partitions){
    double cutSize = 0;
    for(auto& partition : partitions)
        for(auto& node : partition.nodes){
            auto nNodes = G.getNeighbors(node);
            for(auto& nNode:nNodes)
                //check if the neighbour is in the same partition
                if(std::find(partition.nodes.begin(), partition.nodes.end(), nNode) == partition.nodes.end())
                    cutSize += G.getEdgeWeight(node, nNode);
        }
    return cutSize/2;
}
/*
 * Insert to the partition some weight 0 nodes to equal the partition sizes 
 */
void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight){
    for(int i = 0; i < static_cast<int>(max_weight-partition.weight); i++){
        alreadyIns++;
        partition.nodes.push_back(-alreadyIns);
    }
}

/*
 * Divide randomly the nodes into the partitions based on the total number of nodes
 */
std::vector<Partition> makeNodePartion(const Graph& G, const int& numPartitions){
    auto nodes = G.getVertices();
    std::vector<Partition> partitions(numPartitions);
    int i = 0;
    for(auto& node : nodes){
        partitions[i].nodes.push_back(node.first);
        i = (i+1)%numPartitions;
    }
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
    return partitions;
}

/*
 * Divide the nodes in partitions based on the weight
 */
std::vector<Partition> makePartion(Graph& G, const int& numPartitions, int& inserted){
    std::vector<Partition> partitions(numPartitions);
    double totVertexWeight = 0;
    // List of Vertices
    auto vertices = G.getVertices();
    //std::cout << "Number of Vertices :" << vertices.size()<< std::endl;
    /*for(auto vertex: vertices)
        std::cout << vertex.first << "," <<vertex.second << "  ";
    std::cout << std::endl;*/
    // Ordered map to order the node based on their weight
    std::multimap<double, int> nodes;
    for(auto& node: vertices){
        totVertexWeight += node.second;
        nodes.insert(std::pair<double, int>(node.second, node.first));
    }
    // Free first vertices data 
    vertices.clear();
    // Average
    double nodeSize = totVertexWeight/static_cast<double>(numPartitions);
    double excessSize = static_cast<int>(totVertexWeight)%numPartitions;
    std::cout << "Average Node Size: " << nodeSize << std::endl;
    size_t i = 0;
    // divide nodes based on weight
    for(auto it = nodes.rbegin(); it != nodes.rend(); ++it){
        if(partitions[i].weight >= static_cast<int>(nodeSize)){
            size_t j = (i+1)%numPartitions;
            //if average weight reached check next partition until a full cycle or not full one found
            while(partitions[j].weight >= static_cast<int>(nodeSize) && i != j) {
                j = (j+1)%numPartitions;
            }
            i = j;
        }
        partitions[i].nodes.push_back(it ->second);
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
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
    //
    printPartitions(partitions);
    inserted = 0;
    for(auto& partition: partitions){
        std::vector<int> partitionInflatedNodes;
        for(const auto& node : partition.nodes){
            std::vector<int> inflatedNodes = G.inflateVertex(node, inserted);
            partitionInflatedNodes.insert(partitionInflatedNodes.end(), inflatedNodes.begin(), inflatedNodes.end());
        }
        partition.nodes.insert(partition.nodes.begin(), partitionInflatedNodes.begin(), partitionInflatedNodes.end());
    }
    printPartitions(partitions);
    return partitions;
    
}



// std::vector<std::vector<int>> makePartion(const Graph& G, const int& numPartitions){
//     std::vector<std::vector<int>> fPartitions(numPartitions);
//     std::vector<Partition> partitions(numPartitions);
//     double totVertexWeight = 0;
//     auto vertices = G.getVertices();
//     std::cout << "Number of Vertices :" << vertices.size()<< std::endl;
//     /*for(auto vertex: vertices)
//         std::cout << vertex.first << "," <<vertex.second << "  ";
//     std::cout << std::endl;*/
//     std::multimap<double, int> nodes;
//     for(auto& node: vertices){
//         totVertexWeight += node.second;
//         nodes.insert(std::pair<double, int>(node.second, node.first));
//     }
//     vertices.clear();
//     double nodeSize = totVertexWeight/static_cast<double>(numPartitions);
//     double excessSize = static_cast<int>(totVertexWeight)%numPartitions;
//     size_t i = 0;
//     for(auto it = nodes.rbegin(); it != nodes.rend(); ++it){
//         if(partitions[i].weight >= nodeSize){
//             size_t j = (i+1)%numPartitions;
//             while(partitions[j].weight >= nodeSize-0.5 && i != j) {
//                 j = (j+1)%numPartitions;
//             }
//             i = j;
//         }
//         partitions[i].nodes.push_back(it ->second);
//         partitions[i].weight += it->first;
//         i = (i+1)%numPartitions;
//     }
//     /*
//     i = 0; 
//     for(auto p: partitions){
//         i++;
//         std::cout << "Partition " << i << " : weight " << p.weight<< std::endl;
//         for(auto n: p.part)
//             std::cout << n << " ";
//         std::cout << std::endl;
//     }
//     auto tpart = cleanPartition(partitions);
//     for(auto p: tpart){
//         i++;
//         std::cout << "Partition " << i << " :" << std::endl;
//         for(auto n: p)
//             std::cout << n << " ";
//         std::cout << std::endl;
//     }*/
//     std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
//     double max_weight = 0;
//     for(auto& partition: partitions)
//         max_weight = (partition.weight > max_weight) ? partition.weight : max_weight;
//     int inserted = 0;
//     for(auto& partition: partitions)
//         fillPartition(partition, inserted, max_weight);
//     /*
//     for(int j = 0; j < numPartitions; j++) {
//         for(const auto& node: partitions[j].part){
//             std::vector<int> inflatedNodes = G.inflateVertex(node, inserted);
//             fPartitions[j].insert(fPartitions[j].end(), inflatedNodes.begin(), inflatedNodes.end());
//         }
//     }
//     return fPartitions;*/
//     return cleanPartition(partitions);
// }


void fillToSize(std::vector<int>& part, const int goalSize, int& index){
    for(int i = part.size() - goalSize;  i < 0 ; i++){
        index++;
        part.push_back(-index);
    }
}

char balancePartitionsSize(Partition& a, Partition& b, const int& startIndex){
    int i = startIndex;
    if(a.nodes.size() == b.nodes.size()) return 0;
    if(a.nodes.size() > b.nodes.size()){
        fillToSize(b.nodes, a.nodes.size(), i);
    } else if(a.nodes.size() < b.nodes.size()){
        fillToSize(a.nodes, b.nodes.size(), i);
    }
    return 1;
}

void deletePartitionsFiller(Partition& a, Partition& b, const int& fillerIndex){
    a.nodes.erase(std::remove_if(a.nodes.begin(), a.nodes.end(), [fillerIndex](int n) { return n < -fillerIndex; }), a.nodes.end());
    b.nodes.erase(std::remove_if(b.nodes.begin(), b.nodes.end(), [fillerIndex](int n) { return n < -fillerIndex; }), b.nodes.end());
}
void deflate(Graph& G, std::vector<Partition>& partitions){
    G.deflate();
    for(auto& partition: partitions)
        partition.nodes.erase(std::remove_if(partition.nodes.begin(), partition.nodes.end(), [](int n) {return n < 0;}),partition.nodes.end());
}

std::vector<Partition> multipartitionKL(Graph& G, const int& numPartitions) {
    int inflatedNodes;
    auto partitions = makePartion(G, numPartitions, inflatedNodes);
    int fillerNodes = inflatedNodes;
    std::cout << "Number of nodes:" << inflatedNodes << std::endl;
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
                char inequal_size = balancePartitionsSize(partitions[i], partitions[j], inflatedNodes);
                int iterations = KL_Partitioning(G, partitions[i], partitions[j]);
                //int iterations = 0;
                if(inequal_size)
                    deletePartitionsFiller(partitions[i], partitions[j], inflatedNodes);
                if(iterations > 1 && optimized) optimized = 0;
            }
        }
    } while (!optimized);
    std::cout << "Number of Rounds: " << rounds << std::endl;
    deflate(G, partitions);
    for(auto& partition: partitions){
        partition.weight = 0;
        for(auto node: partition.nodes){
            partition.weight += G.getVertexWeight(node);
        }
    }
    return partitions;
}


