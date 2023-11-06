#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include "coarsening.cpp"
#include "uncoarsening.cpp"
#include "Reader.cpp"

using namespace std;

//To be defined in LoadGraphFromMemory
Graph graph;

double getCutSize(Graph& graph, std::vector<std::vector<int>>& partitions){
    std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights = graph.getEdgeWeights();

    double cutSize = 0.0;

    for(const auto& edge : edgeWeights){
        int nodeA = edge.first.first;
        int nodeB = edge.first.second;
        double value = edge.second;

        int partitionA = -1;
        int partitionB = -1;

        for(int i=0; i<partitions.size(); i++){
            if (std::find(partitions[i].begin(), partitions[i].end(), nodeA) != partitions[i].end()) {
                partitionA = i;
            }

            if (std::find(partitions[i].begin(), partitions[i].end(), nodeB) != partitions[i].end()) {
                partitionB = i;
            }
        }

        if(partitionA != partitionB){
            cutSize += value/2;
        }
    }

    return cutSize;
}

std::vector<std::vector<int>> InitialPartitioning(Graph& coarsedGraph, int npartitions, float maxDeviation) {
    // Calculate the total weight of all vertices in the graph
    double total_weight = 0;
    for (const auto& vertexPair : coarsedGraph.getVertices()) {
        total_weight += vertexPair.second;
    }

    // Calculate the target weight for each partition
    double target_weight = total_weight / npartitions;

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

    //Greedy assignment of vertices to partitions
    while(condition){
        int k = -1;

        //if no node can be assigned in one complete cycle, remove the contraint that avoid to add nodes to heavy partitions
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
                    //doing so we try to take into acount max deviation for initial partitioning
                    continue;
                }

                //If the partition is empty or if the node is completely disconnected, assign it
                if(partitions[partitionOrder[i]].empty() || coarsedGraph.getNeighbors(vertex).size() == 0){
                    connected = true;
                    min_weight_partition = partitionOrder[i];
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
            }
        }

        //check if all assigned vector elements are true
        // Check if all values in the bool vector are true
        bool allTrue = std::all_of(assigned.begin(), assigned.end(), [](bool val) { return val; });

        if (allTrue) {
            // If all values are true, set condition to false and exit the loop
            condition = false;
        }
    }

    std::cout << "Start refining initial partitioning" << std::endl;

    int iterations = 0;
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
        for(int i=(partitions.size()-1); i>0; i--){

            if(partitions[partitionOrder[0]].empty()){
                int node = *partitions[partitionOrder[i]].begin();
                //move neighbor to partition
                partitions[partitionOrder[i]].erase(std::remove(partitions[partitionOrder[i]].begin(), partitions[partitionOrder[i]].end(), node), partitions[partitionOrder[i]].end());
                partitions[partitionOrder[0]].push_back(node);
                partition_weights[partitionOrder[i]] -= graph.getVertexWeight(node);
                partition_weights[partitionOrder[0]] += graph.getVertexWeight(node);
                swapped = true;
                break;
            }

            for(int vertex : partitions[partitionOrder[0]]){
                std::vector<int>& neighbors = coarsedGraph.getNeighbors(vertex);

                for(int node : partitions[partitionOrder[i]]){
                    //Check all nodes of the heaviest partition not already checked
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
        iterations++;
    }

    std::cout << "It iterates " << iterations << " times to balance initial partitions" << std::endl;

    return partitions;
}

std::vector<int> findBoudaryVertices(Graph& graph, std::vector<std::vector<int>> initial_partitions) {
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

void WriteOutputToFile(const std::vector<std::vector<int>>& partitions, string outputfile){
    //TO DO
    int partitionIndex = 0;
    for (const auto& partition : partitions) {
        std::cout << "Partition " << partitionIndex << ": ";
        int sum = 0;
        for (const auto& vertex : partition) {
            //std::cout << vertex << " ";
            sum += graph.getVertexWeight(vertex);
        }
        std::cout << std::endl;
        std::cout << "Total weight: " << sum << std::endl;
        std::cout << std::endl;
        partitionIndex++;
    }
}


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile, string outputfile){
    auto start_time = std::chrono::high_resolution_clock::now();
    graph = metisRead(inputfile, nthreads);      //load the graph from file
    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Reading time: " << seconds << " seconds" << std::endl;
    //graph.print();
    start_time = std::chrono::high_resolution_clock::now();
    std::cout << std::endl << std::endl << "COARSENING" << std::endl;
    Graph coarsedGraph = Coarsening(graph, nthreads);        //Coarse the initial graph
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    std::cout << "Coarsening time: " << seconds << " seconds" << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    std::cout << std::endl << "INITIAL PARTITIONING" << std::endl;
    std::vector<std::vector<int>> initial_partitions = InitialPartitioning(coarsedGraph, npartitions, maxdeviation);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    std::cout << "Initial partitioning time: " << seconds << " seconds" << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    std::cout << std::endl << "BOUNDARY VERTICES" << std::endl;
    std::vector<int> boundaryVertices = findBoudaryVertices(coarsedGraph, initial_partitions);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    std::cout << "Boundary vertices time: " << seconds << " seconds" << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    std::cout << std::endl << "UNCOARSENING" << std::endl;
    std::vector<std::vector<int>> uncoarsened_partitions = Uncoarsening(coarsedGraph, initial_partitions, boundaryVertices, nthreads, maxdeviation);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    std::cout << "Uncoarsening time: " << seconds << " seconds" << std::endl;
    WriteOutputToFile(uncoarsened_partitions, outputfile);
    std::cout << "Cut size: " << getCutSize(graph, uncoarsened_partitions) << std::endl;
}


int main() {
    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    //resources\metismodels\x200000y440000m20q20.metis
    const string inputfile = "resources/metismodels/x200000y440000m20q20.metis";
    int nthreads = 10; // Change this to the desired number of threads
    int npartitions = 3; // Change this to the desired number of partitions
    float maxdeviation = 1.05; // Change this to the desired max deviation

    std::string outputfile = "output_partition.txt"; // Change this to the output file name

    // Call the algorithm function
    MultithreadedMETIS(nthreads, npartitions, maxdeviation, inputfile, outputfile);

    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Execution time: " << seconds << " seconds" << std::endl;


    return 0;
}
