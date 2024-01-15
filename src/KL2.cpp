
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
#define NPART 2

// Structure to represent a partition
struct Partition {
    std::vector<int> part;
    double weight;
};

// Function to calculate the cut size of partitions
double calculateCutSize(const Graph& G, const std::vector<std::vector<int>>& partitions){
    double cutSize = 0;
    for(const auto& partition : partitions)
        for(const auto& node : partition){
            auto nNodes = G.getNeighborsKL(node);
            for(const auto& nNode : nNodes)
                if(std::find(partition.begin(), partition.end(), nNode) == partition.end())
                    cutSize += G.getEdgeWeightKL(node, nNode);
        }
    return cutSize / 2;  // Dividing by 2 as each edge is counted twice
}

// Function to fill a partition with dummy nodes to match the maximum weight
void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight){
    DEBUG_STDOUT("filling partition");
    for(int i = 0; i < static_cast<int>(max_weight - partition.weight); i++){
        alreadyIns++;
        partition.part.push_back(-alreadyIns);
    }
}

// Function to print partitions
void printPartitions(const std::vector<std::vector<int>>& partitions, std::unordered_map<int, double> vertices){
    int weight = 0;
    for(int i = 0; i < partitions.size(); i++){
        std::cout << "Partition " << i + 1 << ": " << partitions[i].size() << std::endl;
        for(const auto& node: partitions[i]){
            std::cout << node << " ";
            weight += vertices[node];
        }
        std::cout << " weight: " << weight;
        std::cout << std::endl;
    }
}

// Function to clean partitions by removing dummy nodes
std::vector<std::vector<int>> cleanPartition(const std::vector<Partition>& partitions){
    std::vector<std::vector<int>> cleanPartitions;
    for(const auto& partition : partitions)
        cleanPartitions.push_back(partition.part);
    return cleanPartitions;
}

// Function to create initial node partitions
std::vector<std::vector<int>> makeNodePartion(const Graph& G, const int& numPartitions){
    auto nodes = G.getVertices();
    std::vector<std::vector<int>> partitions(numPartitions);
    int i = 0;
    for(const auto& node : nodes){
        partitions[i].push_back(node.first);
        i = (i + 1) % numPartitions;
    }
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
    return partitions;
}


// Function to create initial node partitions
std::vector<std::vector<int>> makeNodePartion2(const Graph& G, const int& numPartitions){
    auto nodes = G.getVertices();
    auto expandedRange = G.getExpandedRange();
    std::vector<std::vector<int>> partitions(numPartitions);
    int expandedStart = G.getExpandedStart();
    
    // Step 1: Randomly assign each original vertex to partition A or B
    for(const auto& entry : nodes){
        int vertexID = entry.first;
        double weight = entry.second;
        int exStartIndex = expandedRange[vertexID];
        
        int val = rand() % 2 == 0;
        if (vertexID < expandedStart) {
          partitions[val].push_back(entry.first); // put node in random partition
          for(int i; i < weight; i++) {
            partitions[val].push_back(exStartIndex+i); // Step 2: Assign all extra generated vertices to the same partition as their original vertex
          }
        }
    }
    
    
    // // Step 3: Assign about 1/2 of the extra nodes from the largest partition to the partition with fewer nodes
    // int count_partition_a = std::count_if(nodes.begin(), nodes.end(), [](const auto& node){ return node.partition_label == "A"; });
    // int count_partition_b = std::count_if(nodes.begin(), nodes.end(), [](const auto& node){ return node.partition_label == "B"; });
    // std::cout << "Partition A unbalanced weight: " << count_partition_a << std::endl;
    // std::cout << "Partition B unbalanced weight: " << count_partition_b << std::endl;

    // while (count_partition_a > count_partition_b) {
    //     auto extra_node = std::find_if(nodes.begin(), nodes.end(),
    //                                     [](const auto& node){ return node.partition_label == "A" && "_" in std::to_string(node.id); });
    //     if (extra_node != nodes.end()) {
    //         extra_node->partition_label = "B";
    //         count_partition_a--;
    //         count_partition_b++;
    //     }
    // }

    // while (count_partition_b > count_partition_a) {
    //     auto extra_node = std::find_if(nodes.begin(), nodes.end(),
    //                                     [](const auto& node){ return node.partition_label == "B" && "_" in std::to_string(node.id); });
    //     if (extra_node != nodes.end()) {
    //         extra_node->partition_label = "A";
    //         count_partition_a++;
    //         count_partition_b--;
    //     }
    // }

    // std::cout << "Partition A balanced weight: " << count_partition_a << std::endl;
    // std::cout << "Partition B balanced weight: " << count_partition_b << std::endl;

    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
    return partitions;
}

// Function for Kernighan-Lin Multi-level Partitioning
std::vector<std::vector<int>> multipartitionKL(Graph& G, const int& numPartitions) {
    auto partitions = makeNodePartion2(G, numPartitions);
    char optimized = 1;
    int rounds = 0;
    do {
        rounds++;
        cout << "New round" << endl;
        optimized = 1;
        for(int i = 0; i < numPartitions; i++){
            for(int j = i + 1; j < numPartitions; j++){
                cout << "New partitioning " << i << ":" << j << endl;
                printPartitions(partitions, G.getVertices());
                char inequal_size = 0;
                if(partitions[i].size() < partitions[j].size()){
                    partitions[i].push_back(-1);
                    inequal_size = 1;
                } else if (partitions[i].size() > partitions[j].size()){
                    partitions[j].push_back(-1);
                    inequal_size = 1;
                }
                int iterations = KL_Partitioning(G, partitions[i], partitions[j]);
                cout << "Iterations " << iterations << endl;
                if(inequal_size){
                    partitions[i].erase(std::remove(partitions[i].begin(), partitions[i].end(), -1), partitions[i].end());
                    partitions[j].erase(std::remove(partitions[j].begin(), partitions[j].end(), -1), partitions[j].end());
                }
                if(iterations > 1 && optimized) optimized = 0;
            }
        }
    } while (!optimized);
    std::cout << "Number of Rounds: " << rounds << std::endl;
    return partitions;
}

// Main function
int main() {
    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    // Define the file path and number of threads
    const std::string filename = "/content/q2/resources/metismodels/x15y30m20q20.metis";
    // const std::string filename = "/content/q2/resources/metismodels/x1000y2000m20q20.metis";
    const int numThreads = 1;  // Change the number of threads if needed

    // Read the graph from the file
    Graph graph = metisRead(filename, numThreads);

    // Call multi-level KL partitioning
    auto optPartitions = multipartitionKL(graph, NPART);

    // Print the final partitions and some statistics
    printPartitions(optPartitions, graph.getVertices());

    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Execution time: " << seconds << " seconds" << std::endl;
    // Print the final cut size
    std::cout << "Final CutSize: " << calculateCutSize(graph, optPartitions);

    return 0;
}
