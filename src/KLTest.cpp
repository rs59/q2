#include "KLGroups.cpp"
#include"Reader.cpp"

#include <chrono>
#define NPART 7


int main() {
    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();
    //resources\metismodels\x100000y220000m20q20.metis
    const string filename = "resources/metismodels/x100y200m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    //graph.print();
    //std::cout << "Registred Number of Vertices: " << graph.numVertices()<< std::endl;
    
    auto optPartitions = multipartitionKL(graph, NPART);
    printPartitions(optPartitions);
    // Stop the clock
    std::cout << "Final CutSize: " << calculateCutSize(graph, optPartitions);
    return 0;
}


