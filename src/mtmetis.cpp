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


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile, string outputfile){
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
    Graph coarsedGraph = Coarsening(graph, nthreads, npartitions);        //Coarse the initial graph
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    seconds = duration.count() / 1e6;
    DEBUG_STDOUT("Coarsening time: " + std::to_string(seconds) + " seconds");
    start_time = std::chrono::high_resolution_clock::now();

    DEBUG_STDOUT("\nINITIAL PARTITIONING");
    coarsedGraph.setOriginalVertices(coarsedGraph.size());
    std::cout << coarsedGraph.getExpandedStart() << std::endl;
    std::vector<std::vector<int>> initial_partitions = multipartitionKL_blob(coarsedGraph, npartitions);

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
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " nthreads npartitions maxdeviation inputfile outputfile" << std::endl;
        return 1;
    }

    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    int nthreads = std::stoi(argv[1]);
    int npartitions = std::stoi(argv[2]);
    float maxdeviation = std::stof(argv[3]);
    std::string inputfile = argv[4];
    std::string outputfile = argv[5];

    // Call the algorithm function
    MultithreadedMETIS(nthreads, npartitions, maxdeviation, inputfile, outputfile);

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
