#ifdef DEBUG
#include <iostream>
#endif

#include "KLCore.cpp"

// Structure to represent a partition
struct Partition {
    std::vector<int> part;
    double weight;
};

#ifdef DEBUG
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

// Function to print partitions
void printPartitions(const std::vector<std::vector<int>>& partitions, std::unordered_map<int, double> vertices){
    int weight;
    std::cout << "Printing Partitions: " << std::endl;
    for(int i = 0; i < partitions.size(); i++){
        weight = 0;
        std::cout << "Partition " << i + 1 << ": " << partitions[i].size() << std::endl;
        for(const auto& node: partitions[i]){
            std::cout << node << " ";
            weight += vertices[node];
        }
        std::cout << " weight: " << weight;
        std::cout << std::endl;
    }
}
#endif

// Function to fill a partition with dummy nodes to match the maximum weight
void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight){
    DEBUG_STDOUT("filling partition");
    for(int i = 0; i < static_cast<int>(max_weight - partition.weight); i++){
        alreadyIns++;
        partition.part.push_back(-alreadyIns);
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
std::vector<std::vector<int>> makeNodePartition_multiple_round_robin(const Graph& G, const int& numPartitions){
    auto nodes = G.getVertices();
    std::vector<std::vector<int>> partitions(numPartitions);
    int i = 0;
    for(const auto& node : nodes){
        partitions[i].push_back(node.first);
        i = (i + 1) % numPartitions;
    }
#ifdef DEBUG
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
#endif
    return partitions;
}


// Function to create initial node partitions
std::vector<std::vector<int>> makeNodePartition(Graph& G, bool expand){
    int numPartitions = 2;
    srand(time(NULL));
    auto nodes = G.getVertices();
    std::vector<std::vector<int>> partitions(numPartitions);
    int expandedStart = G.getExpandedStart();

    // G.print();
    
    std::unordered_map<int, int> partition_assignments;
    std::unordered_map<int, int> partition_totalweights;

    // Step 1: Randomly assign each original vertex to partition A or B
    for(const auto& entry : nodes){
        int vertexID = entry.first;
        double weight = entry.second;
        DEBUG_STDOUT("id: " + std::to_string(vertexID) + " weight: " + std::to_string(weight) + " exStartIndex: ");

        int val = rand() % 2 == 0;
        if (vertexID < expandedStart + 1)
        {                                           // if it is not an expanded node
            partitions[val].push_back(entry.first); // put node in random partition
            partition_assignments[vertexID] = val;
            partition_totalweights[val] += weight;
        }
    }


    // Step 2: Rebalance partitions in terms of weight: move over 1/2 of nodes from lighter partition
    int lighter_partition;
    int heavier_partition;
    DEBUG_STDOUT("p1 total weights " + std::to_string(partition_totalweights[0]) + "p2 total weights " + std::to_string(partition_totalweights[1]));
    if(partition_totalweights[0] < partition_totalweights[1]) {
        lighter_partition = 0;
    } else {
        lighter_partition = 1;
    }
    heavier_partition = 1-lighter_partition;

    int to_move_minimum = (partition_totalweights[heavier_partition] - partition_totalweights[lighter_partition])/2;
    DEBUG_STDOUT("tmm " + std::to_string(to_move_minimum));
    // Remove nodes from the heavier partition until it is balanced
    int orig_size = partitions[heavier_partition].size();
    for (int i = orig_size - 1; i >= 0; i -= 1)
    {
        int vertexID = partitions[heavier_partition][i];

        int thisWeight = nodes[vertexID];
        to_move_minimum -= thisWeight;
        if (to_move_minimum <= 0)
        {
            break;
        }

        partitions[heavier_partition].pop_back();
        partitions[lighter_partition].push_back(vertexID);
        
        DEBUG_STDOUT("moving " + std::to_string(vertexID) + " from " + std::to_string(heavier_partition) + " to " + std::to_string(lighter_partition) + " with weight " + std::to_string(thisWeight) );
        DEBUG_STDOUT("heavier partition " + std::to_string(partitions[heavier_partition].size()) + " " + std::to_string(heavier_partition) + " and " + std::to_string(partitions[lighter_partition].size()) + " lighter partition " + std::to_string(lighter_partition));
        // std::cout << "moving " << vertexID << " from " << heavier_partition << " to " << lighter_partition << " with weight " << thisWeight << std::endl;
        // std::cout << "heavier partition " << partitions[heavier_partition].size() << " " << heavier_partition << " and " << partitions[lighter_partition].size() << " lighter partition " << lighter_partition << std::endl;

        partition_totalweights[heavier_partition] -= thisWeight;
        partition_totalweights[lighter_partition] += thisWeight;
        DEBUG_STDOUT("lighter p: " + std::to_string(partition_totalweights[lighter_partition]) + "heavier p: " + std::to_string(partition_totalweights[heavier_partition]));
        // std::cout << "lighter p: " << partition_totalweights[lighter_partition] << "heavier p: " << partition_totalweights[heavier_partition] << std::endl;
        partition_assignments[vertexID] = lighter_partition;
#ifdef DEBUG
        std::for_each(partitions[heavier_partition].begin(), partitions[heavier_partition].end(), [](int value)
                    { std::cout << value << " "; });
        std::cout << std::endl;
        std::for_each(partitions[lighter_partition].begin(), partitions[lighter_partition].end(), [](int value)
                    { std::cout << value << " "; });
        std::cout << std::endl;
#endif
    }

    // We can try to expand the nodes but it seems to be very slow and not worth it
    if(expand) {
        G.setOriginalVertices(G.size());
        DEBUG_STDOUT(std::to_string(G.size()) + " VERTICES IN GRAPH");
        // std::cout << G.size() << " VERTICES IN GRAPH" << std::endl;

        G.expandNodes();
        auto expandedRange = G.getExpandedRange();
        

        // Assign expanded nodes t
        for(const auto& entry : nodes){
            int vertexID = entry.first;
            double weight = entry.second;
            int exStartIndex = expandedRange[vertexID];
            
            if (vertexID < expandedStart+1) { // if it is not an expanded node
                int val = partition_assignments[vertexID];
                
                for(int i=0; i < weight; i++) {  // if it is an expanded node
                    // cout << "val : " << val << "currInd : " << exStartIndex+i << " expandedRange : " << expandedRange[exStartIndex+i] << std::endl;
                    partitions[val].push_back(exStartIndex+i); // Step 2: Assign all extra generated vertices to the same partition as their original vertex
                }
            
            }
        }
        
        
        // De-expand nodes

        G.collapseNodes();

    }

    // std::cout << "Partition A balanced weight: " << count_partition_a << std::endl;
    // std::cout << "Partition B balanced weight: " << count_partition_b << std::endl;
#ifdef DEBUG
    std::cout << "Initial CutSize: " << calculateCutSize(G, partitions) << std::endl;
#endif
    return partitions;
}


std::pair<int, int> dividePartitionCount(const int &n)
{
    int part1;
    int part2;
    if (n % 2 == 1)
    {
        part1 = n / 2;
        part2 = n / 2 + 1;
    }
    else
    {
        part1 = part2 = n / 2;
    }
    return std::make_pair(part1, part2);
}

// Function for Kernighan-Lin Multi-level Partitioning
std::vector<std::vector<int>> multipartitionKL_blob(Graph &G, int numPartitions)
{

    int gOriginalSize = G.size();
    G.setOriginalVertices(G.size());
    // std::pair<int, int> partitionsToMake = dividePartitionCount(numPartitions);
    std::vector<std::vector<int>> partitions = makeNodePartition(G, false);
    numPartitions--; // one less split to complete
    if (numPartitions == 1)
    { // if there is only one partition left to split, then quit
        return partitions;
    }
    std::vector<Graph> graphs;
    std::pair<Graph, Graph> graphPair = G.splitGraph(partitions[0], partitions[1]);
    graphs.push_back(graphPair.first);
    graphs.push_back(graphPair.second);

    // Now partitions and graphs match and are in vectors
#ifdef DEBUG
    printPartitions(partitions, G.getVertices());
#endif
    while (true) // Divide all of the partitions in the last round (2,4,8,16,etc.)
    {

        auto currentGraphs = graphs;
        int orig_size = partitions.size();
        for (auto tempG : currentGraphs)
        {
            DEBUG_STDOUT("number of graphs: " + std::to_string(graphs.size()));
            // std::cout << "number of graphs: " << graphs.size() << std::endl;
            // tempG.print();
            // splitGraph for each partition
            // std::cout << "tempG size " << tempG.size() << std::endl;

            tempG.setOriginalVertices(gOriginalSize);
            auto tempPartitions = makeNodePartition(tempG, false);
            KL_Partitioning(tempG, tempPartitions[0], tempPartitions[1]);

            // std::for_each(tempPartitions[0].begin(), tempPartitions[0].end(), [](int value)
            //               { std::cout << value << " "; });
            // std::cout << std::endl;
            // std::for_each(tempPartitions[1].begin(), tempPartitions[1].end(), [](int value)
            //               { std::cout << value << " "; });
            // std::cout << std::endl;

            std::pair<Graph, Graph> tempGraphs = tempG.splitGraph(tempPartitions[0], tempPartitions[1]);
            // printPartitions(tempPartitions, tempGraphs.first.getVertices());

            partitions.erase(partitions.begin());
            partitions.push_back(tempPartitions[0]);
            partitions.push_back(tempPartitions[1]);

            graphs.erase(graphs.begin());
            graphs.push_back(tempGraphs.first);
            graphs.push_back(tempGraphs.second);

            // std::cout << "NEW number of graphs: " << graphs.size() << std::endl;

            numPartitions--; // one less split to complete
            if (numPartitions == 1)
            { // if there is only one partition left to split, then quit
                break;
            }
        }
        if (numPartitions == 1) {
            break;
        }
        
    }
    return partitions;
}

// Function for Kernighan-Lin Multi-level Partitioning
std::vector<std::vector<int>> multipartitionKL_round_robin(Graph& G, const int& numPartitions) {
    auto partitions = makeNodePartition_multiple_round_robin(G, numPartitions);
    char optimized = 1;
    int rounds = 0;
    do {
        rounds++;
        std::cout << "New round" << std::endl;
        optimized = 1;
        for(int i = 0; i < numPartitions; i++){
            for(int j = i + 1; j < numPartitions; j++){
                std::cout << "New partitioning " << i << ":" << j << std::endl;
#ifdef DEBUG                
                printPartitions(partitions, G.getVertices());
#endif
                char inequal_size = 0;
                // Fill with weight 0 nodes to equalize partitions size
                if(partitions[i].size() < partitions[j].size()){
                    partitions[i].push_back(-1);
                    inequal_size = 1;
                } else if (partitions[i].size() > partitions[j].size()){
                    partitions[j].push_back(-1);
                    inequal_size = 1;
                }
                int iterations = KL_Partitioning(G, partitions[i], partitions[j]);
                DEBUG_STDOUT("Iterations " + std::to_string(iterations));
                // std::cout << "Iterations " << iterations << std::endl;
                //Delete weight 0 nodes
                if(inequal_size){
                    partitions[i].erase(std::remove(partitions[i].begin(), partitions[i].end(), -1), partitions[i].end());
                    partitions[j].erase(std::remove(partitions[j].begin(), partitions[j].end(), -1), partitions[j].end());
                }
                if(iterations > 1 && optimized) optimized = 0;
            }
        }
    } while (!optimized);
    DEBUG_STDOUT("Number of Rounds: " +  std::to_string(rounds));
#ifdef DEBUG
    printPartitions(partitions, G.getVertices());
#endif
    // std::cout << "Number of Rounds: " << rounds << std::endl;
    return partitions;
}