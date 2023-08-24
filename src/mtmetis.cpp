#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <thread>
#include <mutex>
#include "coarsening.cpp"
#include "uncoarsening.cpp"

using namespace std;

int NUM_ITERATIONS = 100; //TO BE DEFINED

//To be defined in LoadGraphFromMemory
Graph graph;
Graph coarsedGraph;

//Data structures and variables used in parallel refinement step
std::unordered_map<std::pair<int, int>, double, HashPair> edge_weights;
std::vector<double> vertex_gains;
int num_vertices;
int num_partitions;
double total_weight;
std::unordered_set<int> locked_vertices;  //vertices that won't be moved in this iteration

void LoadGraphFromMemory(string inputfile){
    //TO-DO
}

std::vector<std::vector<int>> InitialPartitioning(int npartitions) {
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

    // Greedy assignment of vertices to partitions
    for (int vertex : sorted_vertices) {
        // Find the partition with the minimum current weight
        int min_weight_partition = 0;
        double min_weight = partition_weights[0];
        for (int i = 1; i < npartitions; ++i) {
            if (partition_weights[i] < min_weight) {
                min_weight_partition = i;
                min_weight = partition_weights[i];
            }
        }

        // Add the vertex to the partition
        partitions[min_weight_partition].push_back(vertex);
        partition_weights[min_weight_partition] += coarsedGraph.getVertexWeight(vertex);
    }

    return partitions;
}

// Function to perform a single step of the refinement process using the Kernighan-Lin algorithm
// to improve the partition quality.
void RefinementStep(std::mutex& mtx, std::vector<std::vector<int>>& partition, const int start, const int end, const int n_it) {
    int vertex_partition;

    // Calculate the initial gain for each vertex in its current partition
    for (int vertex = start; vertex < end; vertex++) {
        vertex_partition = -1;
        for (int i = 0; i < num_partitions; i++) {
            if (std::find(partition[i].begin(), partition[i].end(), vertex) != partition[i].end()) {
                vertex_partition = i;
                break;
            }
        }

        for (int neighbor = 0; neighbor < num_vertices; neighbor++) {
            if (neighbor != vertex) {
                auto edge_weight_it = edge_weights.find({vertex, neighbor});
                double weight = (edge_weight_it != edge_weights.end()) ? edge_weight_it->second : 0.0;

                if (std::find(partition[vertex_partition].begin(), partition[vertex_partition].end(), neighbor) !=
                    partition[vertex_partition].end()) {
                    vertex_gains[vertex] -= weight;
                } else {
                    vertex_gains[vertex] += weight;
                }
            }
        }
    }

    // Perform swaps to improve the partition quality
    for (int iter = 0; iter < n_it; iter++) { // A heuristic to limit the number of swaps

        // Find the best pair of vertices to swap
        double best_gain = 0.0;
        int vertex_to_move = -1;
        int vertex_to_stay = -1;

        for (int vertex = start; vertex < end; vertex++) {
            if (locked_vertices.find(vertex) == locked_vertices.end() && (vertex_to_move == -1 || vertex_gains[vertex] > best_gain)) {
                vertex_to_move = vertex;
                best_gain = vertex_gains[vertex];
            }
        }

        locked_vertices.insert(vertex_to_move);

        int neighbor_partition;

        vector<int> neighbors = graph.getNeighbors(vertex_to_move);
        for(int i=0; i < neighbors.size(); i++){
            int neighbor = neighbors[i];
            neighbor_partition = -1;
            for (int i = 0; i < num_partitions; i++) {
                if (std::find(partition[i].begin(), partition[i].end(), neighbor) != partition[i].end()) {
                    neighbor_partition = i;
                    break;
                }
            }

            double gain = vertex_gains[vertex_to_move];
            gain -= 2.0 * graph.getEdgeWeight(vertex_to_move, neighbor);

            if (gain > best_gain) {
                vertex_to_stay = neighbor;
                best_gain = gain;
            }
        }

        locked_vertices.insert(vertex_to_stay);

        //Adjusted until here

        // Swap the two selected vertices between their partitions
        for (int i = 0; i < num_partitions; i++) {
            auto it = std::find(partition[i].begin(), partition[i].end(), vertex_to_move);
            if (it != partition[i].end()) {
                partition[i].erase(it);
                break;
            }
        }

        for (int i = 0; i < num_partitions; i++) {
            auto it = std::find(partition[i].begin(), partition[i].end(), vertex_to_stay);
            if (it != partition[i].end()) {
                partition[i].erase(it);
                break;
            }
        }

        partition[neighbor_partition - start].push_back(vertex_to_move);
        partition[vertex_partition - start].push_back(vertex_to_stay);

        // Update the gains for the affected vertices
        for (int neighbor = 0; neighbor < num_vertices; neighbor++) {
            if (neighbor != vertex_to_move) {
                if (std::find(partition[neighbor_partition].begin(), partition[neighbor_partition].end(), neighbor) != partition[neighbor_partition].end()) {
                    vertex_gains[neighbor] -= 2.0 * graph.getEdgeWeight(vertex_to_move, neighbor);
                }
                else {
                    vertex_gains[neighbor] += 2.0 * graph.getEdgeWeight(vertex_to_move, neighbor);
                }
            }
        }
    }
}


// Function to perform the refinement process in a multithreaded manner.
void MultithreadedRefinement(int nthreads, std::vector<std::vector<int>>& initial_partitions) {
    int num_iterations = NUM_ITERATIONS;

    // Calculate the number of vertices in each chunk for each thread.
    int chunk_size = graph.numVertices() / nthreads;

    //initialize the mutex
    std::mutex mtx;

    //Initialization of common variables and data structures for threads refinement step
    edge_weights = graph.getEdgeWeights();

    // Get the number of vertices and partitions
    int num_vertices = graph.numVertices();
    int num_partitions = initial_partitions.size();

    // Calculate the target weight for each partition
    double total_weight = 0.0;
    for (const auto& entry : edge_weights) {
        total_weight += entry.second;
    }
    double target_weight = total_weight / num_partitions;

    // Initialize the partition gain for each vertex
    vertex_gains = std::vector<double>(num_vertices, 0.0);

    // Create and start the threads.
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        int start_vertex = i * chunk_size;
        int end_vertex = (i + 1) * chunk_size;
        if (i == nthreads - 1) {
            end_vertex = graph.numVertices();
        }

        threads.emplace_back([i, start_vertex, end_vertex, num_iterations, &initial_partitions, &mtx]() {
            RefinementStep(std::ref(mtx), initial_partitions, start_vertex, end_vertex, num_iterations);
        });
    }

    // Wait for all threads to finish.
    for (auto& thread : threads) {
        thread.join();
    }
}

void WriteOutputToFile(const std::vector<std::vector<int>>& partitions, string outputfile){
    //TO DO
    int partitionIndex = 0;
    for (const auto& partition : partitions) {
        std::cout << "Partition " << partitionIndex << ": ";
        int sum = 0;
        for (const auto& vertex : partition) {
            std::cout << vertex << " ";
            sum += graph.getVertexWeight(vertex);
        }
        std::cout << std::endl;
        std::cout << "Total weight: " << sum << std::endl;
        std::cout << std::endl;
        partitionIndex++;
    }
}


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile, string outputfile){
    LoadGraphFromMemory(inputfile);      //load the graph from file

    coarsedGraph = Coarsening(graph);        //Coarse the initial graph

    std::vector<std::vector<int>> initial_partitions = InitialPartitioning(npartitions);

    std::vector<std::vector<int>> uncoarsened_partitions = UncoarsePartitions(coarsedGraph, initial_partitions);

    cout << "Write partitions before refinement:" << endl;

    WriteOutputToFile(uncoarsened_partitions, outputfile);

    MultithreadedRefinement(nthreads, uncoarsened_partitions);

    Graph restoredGraph = Uncoarsening(coarsedGraph);

    cout << "Write partitions after refinement:" << endl;

    WriteOutputToFile(uncoarsened_partitions, outputfile);

}


int main() {
    // Add vertices and edges to the global graph instance
    graph.addVertex(0, 10.0);
    graph.addVertex(1, 25.0);
    graph.addVertex(2, 31.0);
    graph.addVertex(3, 24.0);
    graph.addVertex(4, 65.0);
    graph.addVertex(5, 26.0);
    graph.addVertex(6, 77.0);
    graph.addVertex(7, 18.0);
    graph.addVertex(8, 99.0);
    graph.addVertex(9, 43.0);
    graph.addVertex(10, 111.0);
    graph.addVertex(11, 112.0);
    graph.addVertex(12, 54.0);
    graph.addVertex(13, 71.0);
    graph.addVertex(14, 15.0);
    graph.addVertex(15, 1.0);
    graph.addVertex(16, 112.0);
    graph.addVertex(17, 23.0);
    graph.addVertex(18, 15.0);
    graph.addVertex(19, 18.0);
    graph.addVertex(20, 4.0);
    graph.addVertex(21, 65.0);
    graph.addVertex(22, 23.0);
    graph.addVertex(23, 81.0);
    graph.addVertex(24, 99.0);
    graph.addVertex(25, 5.0);
    graph.addVertex(26, 12.0);
    graph.addVertex(27, 11.0);
    graph.addVertex(28, 33.0);
    graph.addVertex(29, 45.0);

    graph.addEdge(0, 1, 0.5);
    graph.addEdge(0, 2, 0.3);
    graph.addEdge(0, 3, 0.6);
    graph.addEdge(0, 4, 0.7);
    graph.addEdge(1, 5, 0.9);
    graph.addEdge(1, 6, 0.8);
    graph.addEdge(2, 7, 0.2);
    graph.addEdge(2, 8, 0.4);
    graph.addEdge(3, 9, 0.1);
    graph.addEdge(3, 10, 0.3);
    graph.addEdge(4, 11, 0.6);
    graph.addEdge(4, 12, 0.8);
    graph.addEdge(5, 13, 0.5);
    graph.addEdge(5, 14, 0.7);
    graph.addEdge(6, 15, 0.2);
    graph.addEdge(6, 16, 0.4);
    graph.addEdge(7, 17, 0.9);
    graph.addEdge(7, 18, 0.1);
    graph.addEdge(8, 19, 0.3);
    graph.addEdge(8, 20, 0.6);
    graph.addEdge(9, 21, 0.8);
    graph.addEdge(9, 22, 0.5);
    graph.addEdge(10, 23, 0.7);
    graph.addEdge(10, 24, 0.2);
    graph.addEdge(11, 25, 0.4);
    graph.addEdge(11, 26, 0.9);
    graph.addEdge(12, 27, 0.1);
    graph.addEdge(12, 28, 0.3);
    graph.addEdge(13, 29, 0.6);
    graph.addEdge(13, 0, 0.8);
    graph.addEdge(14, 1, 0.5);
    graph.addEdge(15, 2, 0.7);
    graph.addEdge(16, 3, 0.4);
    graph.addEdge(17, 4, 0.2);
    graph.addEdge(18, 5, 0.9);
    graph.addEdge(19, 6, 0.1);
    graph.addEdge(20, 7, 0.3);
    graph.addEdge(21, 8, 0.6);
    graph.addEdge(22, 9, 0.8);
    graph.addEdge(23, 10, 0.5);
    graph.addEdge(24, 11, 0.7);
    graph.addEdge(25, 12, 0.4);
    graph.addEdge(26, 13, 0.2);
    graph.addEdge(27, 14, 0.9);
    graph.addEdge(28, 15, 0.1);
    graph.addEdge(29, 16, 0.3);

    // Set the number of threads and partitions
    int nthreads = 1; // Change this to the desired number of threads
    int npartitions = 3; // Change this to the desired number of partitions
    float maxdeviation = 1.1; // Change this to the desired max deviation
    std::string inputfile = "input_graph.txt"; // Change this to the input file name
    std::string outputfile = "output_partition.txt"; // Change this to the output file name

    // Call the algorithm function
    MultithreadedMETIS(nthreads, npartitions, maxdeviation, inputfile, outputfile);

    cout << "THE END" << endl;

    return 0;
}
