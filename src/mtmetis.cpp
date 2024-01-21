#define DEBUG 1

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include "coarsening.cpp"
#include "uncoarsening.cpp"
#include "ReaderWriter2.cpp"
#include "KL2.cpp"

using namespace std;

//To be defined in LoadGraphFromMemory
Graph graph;

std::vector<std::vector<int>> InitialPartitioning(Graph& coarsedGraph, int npartitions, float maxDeviation) {
    // Calculate the total weight of all vertices in the graph
    double total_weight = 0;
    for (const auto& vertexPair : coarsedGraph.getVertices()) {
        total_weight += vertexPair.second;
    }

    // Calculate the target weight for each partition
    double target_weight = total_weight / npartitions;
    DEBUG_STDOUT("Target weight: "+std::to_string(target_weight));

    // Initialize the partitions and their current weight
    std::vector<std::vector<int>> partitions(npartitions);
    std::vector<double> partition_weights(npartitions, 0);

    // Create an index vector to sort vertices based on their weights
    std::vector<int> sorted_vertices;
    for (const auto& vertexPair : coarsedGraph.getVertices()) {
        sorted_vertices.push_back(vertexPair.first);
    }

    // Sort vertices in descending order of their weights
    std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](int a, int b) {
        return coarsedGraph.getVertexWeight(a) > coarsedGraph.getVertexWeight(b);
    });

    //this will be used to check that all vertices has been assigned to a partition
    std::vector<bool> assigned(sorted_vertices.size(), false);
    bool condition = true;
    bool madeAssignemt = true;  //check if the previous turn we made an assignment, if not remove a constraint

    int assignedVerticesCtr = 0;

    bool noNodeConnected = false;
    //Greedy assignment of vertices to partitions
    while(condition){
        int k = -1;

        //if no node can be assigned in one complete cycle, remove the constraint that avoid to add nodes to heavy partitions
        bool useConstraint = madeAssignemt;
        madeAssignemt = false;

        for(int vertex : sorted_vertices){
            k++;   //this will be used ony to update assigned vector

            if(assigned[k] == true){
                //node is already assigned
                continue;
            }
            //we try to find a partition connected to the node
            bool connected = false;
            int min_weight_partition = 0;

            std::vector<int> partitionOrder(partition_weights.size());

            // Initialize partitionOrder to contain indices from 0 to n-1
            for (int i = 0; i < partitionOrder.size(); ++i) {
                partitionOrder[i] = i;
            }

            // Sort partitionOrder based on partitionWeights in ascending order
            std::sort(partitionOrder.begin(), partitionOrder.end(), [&](int a, int b) {
                return partition_weights[a] < partition_weights[b];
            });

            //check if node is connected starting checking the partition with lower weight, and proceding in increasing weight order
            for(int i=0; i < npartitions; i++){
                if(partition_weights[partitionOrder[i]] > maxDeviation*target_weight && useConstraint){
                    //doing so we try to take into account max deviation for initial partitioning
                    continue;
                }

                //If the partition is empty or if the node is completely disconnected or if none of the remaining nodes are connected to any partition, assign it
                if(partitions[partitionOrder[i]].empty() || coarsedGraph.getNeighbors(vertex).size() == 0 || noNodeConnected){
                    connected = true;
                    min_weight_partition = partitionOrder[i];
                    noNodeConnected = false;
                    break;
                }
                for (int node : partitions[partitionOrder[i]]) {
                    if (coarsedGraph.containsEdge(vertex, node)) {
                        connected = true;
                        min_weight_partition = partitionOrder[i];
                        break;
                    }
                }
                if(connected == true){
                    break;
                }
            }
            if(connected){
                assigned[k] = true;
                partitions[min_weight_partition].push_back(vertex);
                partition_weights[min_weight_partition] += coarsedGraph.getVertexWeight(vertex);
                madeAssignemt = true;
                assignedVerticesCtr++;
            }
        }

        if(madeAssignemt == false && useConstraint == false){
            //All nodes in all partitions are disconnected from the remaining nodes
            noNodeConnected = true;
        }

        //check if all assigned vector elements are true
        // Check if all values in the bool vector are true
        bool allTrue = std::all_of(assigned.begin(), assigned.end(), [](bool val) { return val; });

        if (allTrue) {
            // If all values are true, set condition to false and exit the loop
            condition = false;
        } else {
            DEBUG_STDOUT("Condition not met");
            std::cout << "Assigned nodes " << assignedVerticesCtr << " over " << coarsedGraph.getVertices().size() << std::endl;
        }
    }

    int iterations = 0;
    bool constraint = false;
    int constraintCtr = 0;
    //Perform movements until is balanced, or until it iterates one time for vertex
    while(true && iterations < coarsedGraph.getVertices().size()){
        bool balanced = true;
        for(double pw : partition_weights){
            if(pw > maxDeviation*target_weight || pw < (1 - (maxDeviation - 1))*target_weight){
                balanced = false;
                break;
            }
        }
        if(balanced){
            break;
        }

        bool swapped = false;

        std::vector<int> partitionOrder(partition_weights.size());

        // Initialize partitionOrder to contain indices from 0 to n-1
        for (int i = 0; i < partitionOrder.size(); ++i) {
            partitionOrder[i] = i;
        }

        // Sort partitionOrder based on partitionWeights in descending order
        std::sort(partitionOrder.begin(), partitionOrder.end(), [&](int a, int b) {
            return partition_weights[a] < partition_weights[b];
        });

        //for each vertex in lightest partition
        for(int i=(partitions.size()-1); i>=(partitions.size()-1); i--){
            int node;
            if(partitions[partitionOrder[0]].empty()){
                if(!partitions[partitionOrder[i]].empty()){
                    node = *partitions[partitionOrder[i]].begin();
                }else{
                    continue;
                }

                //move neighbor to partition
                partitions[partitionOrder[i]].erase(std::remove(partitions[partitionOrder[i]].begin(), partitions[partitionOrder[i]].end(), node), partitions[partitionOrder[i]].end());
                partitions[partitionOrder[0]].push_back(node);
                partition_weights[partitionOrder[i]] -= coarsedGraph.getVertexWeight(node);
                partition_weights[partitionOrder[0]] += coarsedGraph.getVertexWeight(node);
                swapped = true;
                break;
            }

            for(int vertex : partitions[partitionOrder[0]]){
                std::vector<int>& neighbors = coarsedGraph.getNeighbors(vertex);

                for(int node : partitions[partitionOrder[i]]){
                    //Check all nodes of the heaviest partition not already checked
                    if(coarsedGraph.getNeighbors(node).size() == 0 || constraint){
                        //Evaluated node in heaviest partition is disconnected, swap it
                        partitions[partitionOrder[i]].erase(std::remove(partitions[partitionOrder[i]].begin(), partitions[partitionOrder[i]].end(), node), partitions[partitionOrder[i]].end());
                        partitions[partitionOrder[0]].push_back(node);
                        partition_weights[partitionOrder[i]] -= coarsedGraph.getVertexWeight(node);
                        partition_weights[partitionOrder[0]] += coarsedGraph.getVertexWeight(node);
                        swapped = true;
                        constraint = false;
                        break;
                    }

                    for(int neighbor : neighbors){
                        if(neighbor == node){
                            //move neighbor to partition
                            partitions[partitionOrder[i]].erase(std::remove(partitions[partitionOrder[i]].begin(), partitions[partitionOrder[i]].end(), neighbor), partitions[partitionOrder[i]].end());
                            partitions[partitionOrder[0]].push_back(node);
                            partition_weights[partitionOrder[i]] -= coarsedGraph.getVertexWeight(neighbor);
                            partition_weights[partitionOrder[0]] += coarsedGraph.getVertexWeight(neighbor);
                            swapped = true;
                            break;
                        }
                    }
                    if(swapped){
                        break;
                    }
                }
                if(swapped){
                    break;
                }
            }
            if(swapped){
                break;
            }
        }
        if(!swapped){
            constraint = true;
            constraintCtr++;
        }
        iterations++;
    }

    return partitions;
}

std::vector<int> findBoundaryVertices(Graph& graph, std::vector<std::vector<int>> initial_partitions) {
    std::vector<int> boundaryVertices;
    std::unordered_map<int, double> vertices = graph.getVertices();
    std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights = graph.getEdgeWeights();

    // Iterate through each vertex
    for (const auto& vertexEntry : vertices) {
        int vertexID = vertexEntry.first;
        const std::vector<int>& neighbors = graph.getNeighbors(vertexID);

        // Get the partition of the current vertex
        int partitionVertex = -1; // Initialize to an invalid partition

        // Find the partition that contains the current vertex
        for (int i = 0; i < initial_partitions.size(); ++i) {
            if (std::find(initial_partitions[i].begin(), initial_partitions[i].end(), vertexID) != initial_partitions[i].end()) {
                partitionVertex = i;
                break;
            }
        }

        // Initialize a flag to check if the vertex is a boundary vertex
        bool isBoundary = false;

        if(neighbors.size() == 0){
            //isolated node, add to boundaries
            isBoundary = true;
        }

        // Iterate through the neighbors of the current vertex
        for (int neighbor : neighbors) {
            // Find the partition that contains the neighbor
            int partitionNeighbor = -1; // Initialize to an invalid partition
            for (int i = 0; i < initial_partitions.size(); ++i) {
                if (std::find(initial_partitions[i].begin(), initial_partitions[i].end(), neighbor) != initial_partitions[i].end()) {
                    partitionNeighbor = i;
                    break;
                }
            }

            // Check if the neighbor is in a different partition
            if (partitionVertex != partitionNeighbor) {
                isBoundary = true;
                break;  // No need to check further if it's already a boundary vertex
            }
        }

        // If the vertex is a boundary vertex, add it to the list
        if (isBoundary) {
            boundaryVertices.push_back(vertexID);
        }
    }

    return boundaryVertices;
}

void PrintDetails(const std::vector<std::vector<int>>& partitions, string outputfile){
    //TO DO
    int partitionIndex = 0;
    for (const auto& partition : partitions) {
        DEBUG_STDOUT("Partition " + std::to_string(partitionIndex) + ": ");
        int sum = 0;
        for (const auto& vertex : partition) {
            sum += graph.getVertexWeight(vertex);
        }
        DEBUG_STDOUT("Total weight: " + std::to_string(sum) + "\n");
        partitionIndex++;
    }
}


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile, string outputfile, int coarsestSize, string algorithm){
    auto start_time = std::chrono::high_resolution_clock::now();
    DEBUG_STDOUT("Entering metisRead");
    graph = metisRead(inputfile, nthreads);      //load the graph from file

    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    DEBUG_STDOUT("Reading time: " + std::to_string(seconds) + " seconds");

    start_time = std::chrono::high_resolution_clock::now();
    DEBUG_STDOUT("\nCOARSENING");
    Graph coarsedGraph = Coarsening(graph, nthreads, npartitions, coarsestSize);        //Coarse the initial graph
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    DEBUG_STDOUT("Coarsening time: " + std::to_string(seconds) + " seconds");
    start_time = std::chrono::high_resolution_clock::now();

    DEBUG_STDOUT("\nINITIAL PARTITIONING");
    std::vector<std::vector<int>> initial_partitions;
    if(algorithm == "-g"){
        initial_partitions = InitialPartitioning(coarsedGraph, npartitions, maxdeviation);
    }else{
        coarsedGraph.setOriginalVertices(coarsedGraph.size());
        std::cout << coarsedGraph.getExpandedStart() << std::endl;
        initial_partitions = multipartitionKL_blob(coarsedGraph, npartitions);
    }


    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    DEBUG_STDOUT("Initial partitioning time: " + std::to_string(seconds) + " seconds");
    start_time = std::chrono::high_resolution_clock::now();
    DEBUG_STDOUT("\nBOUNDARY VERTICES");
    std::vector<int> boundaryVertices = findBoundaryVertices(coarsedGraph, initial_partitions);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    DEBUG_STDOUT("Boundary vertices time: " + std::to_string(seconds) + " seconds");
    start_time = std::chrono::high_resolution_clock::now();
    DEBUG_STDOUT("\nUNCOARSENING");
    std::vector<std::vector<int>> uncoarsened_partitions = Uncoarsening(coarsedGraph, initial_partitions, boundaryVertices, nthreads, maxdeviation);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    DEBUG_STDOUT("Uncoarsening time: " + std::to_string(seconds) + " seconds");
    writeToFile(uncoarsened_partitions, outputfile);
    PrintDetails(uncoarsened_partitions, outputfile);
}



int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " nthreads npartitions maxdeviation inputfile outputfile coarsestGraphSize partitioningAlg: -g for greedy or -kl for KL" << std::endl;
        return 1;
    }

    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    int nthreads = std::stoi(argv[1]);
    int npartitions = std::stoi(argv[2]);
    float maxdeviation = std::stof(argv[3]);
    std::string inputfile = argv[4];
    std::string outputfile = argv[5];
    int coarsestSize = std::stoi(argv[6]);
    std::string algorithm = argv[7];

    if(algorithm != "-g" && algorithm != "-kl"){
        std::cerr << "invalid partitioningAlg or partitioningAlg not selected" << std::endl;
        exit(1);
    }

    // Call the algorithm function
    MultithreadedMETIS(nthreads, npartitions, maxdeviation, inputfile, outputfile, coarsestSize, algorithm);

    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    DEBUG_STDOUT("Execution time: " + std::to_string(seconds) + " seconds");


    return 0;
}
