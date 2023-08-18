#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <thread>
#include "coarsening.cpp"
#include "uncoarsening.cpp"

using namespace std;

int NUM_ITERATIONS = 100; //TO BE DEFINED

//To be defined in LoadGraphFromMemory
Graph graph;
Graph coarsedGraph;

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
    std::sort(sorted_vertices.begin(), sorted_vertices.end(), [](int a, int b) {
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
void RefinementStep(std::vector<std::vector<int>>& partition, const int i1, const int i2, const int i3, std::vector<std::vector<int>>& thread_partitions) {
    std::unordered_map<std::pair<int, int>, double, HashPair> edge_weights = coarsedGraph.getEdgeWeights();

    // Get the number of vertices and partitions
    int num_vertices = coarsedGraph.numVertices();
    int num_partitions = partition.size();

    // Calculate the target weight for each partition
    double total_weight = 0.0;
    for (const auto& entry : edge_weights) {
        total_weight += entry.second;
    }
    double target_weight = total_weight / num_partitions;

    // Initialize the partition gain for each vertex
    std::vector<double> vertex_gains(num_vertices, 0.0);

    int vertex_partition;

    // Calculate the initial gain for each vertex in its current partition
    for (int vertex = i1; vertex < i2; ++vertex) {
        vertex_partition = -1;
        for (int i = 0; i < num_partitions; ++i) {
            if (std::find(partition[i].begin(), partition[i].end(), vertex) != partition[i].end()) {
                vertex_partition = i;
                break;
            }
        }

        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex) {
                auto edge_weight_it = edge_weights.find({vertex, neighbor});
                double weight = (edge_weight_it != edge_weights.end()) ? edge_weight_it->second : 0.0;

                if (vertex_partition == -1) {
                    vertex_gains[vertex] += weight;
                }
                else if (std::find(partition[vertex_partition].begin(), partition[vertex_partition].end(), neighbor) != partition[vertex_partition].end()) {
                    vertex_gains[vertex] -= weight;
                }
                else {
                    vertex_gains[vertex] += weight;
                }
            }
        }
    }

    // Initialize the set of locked vertices (vertices that won't be moved in this iteration)
    std::unordered_set<int> locked_vertices;

    // Perform swaps to improve the partition quality
    for (int iter = 0; iter < i3; ++iter) { // A heuristic to limit the number of swaps

        // Find the best pair of vertices to swap
        double best_gain = 0.0;
        int vertex_to_move = -1;
        int vertex_to_stay = -1;

        for (int vertex = i1; vertex < i2; ++vertex) {
            if (locked_vertices.find(vertex) == locked_vertices.end() && (vertex_to_move == -1 || vertex_gains[vertex] > best_gain)) {
                vertex_to_move = vertex;
                best_gain = vertex_gains[vertex];
            }
        }

        locked_vertices.insert(vertex_to_move);

        int neighbor_partition;

        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex_to_move) {
                neighbor_partition = -1;
                for (int i = 0; i < num_partitions; ++i) {
                    if (std::find(partition[i].begin(), partition[i].end(), neighbor) != partition[i].end()) {
                        neighbor_partition = i;
                        break;
                    }
                }

                double gain = vertex_gains[vertex_to_move];
                if (neighbor_partition != -1) {
                    gain -= 2.0 * coarsedGraph.getEdgeWeight(vertex_to_move, neighbor);
                }

                if (gain > best_gain) {
                    vertex_to_stay = neighbor;
                    best_gain = gain;
                }
            }
        }

        locked_vertices.insert(vertex_to_stay);

        // Swap the two selected vertices between their partitions
        for (int i = 0; i < num_partitions; ++i) {
            auto it = std::find(partition[i].begin(), partition[i].end(), vertex_to_move);
            if (it != partition[i].end()) {
                partition[i].erase(it);
                break;
            }
        }

        for (int i = 0; i < num_partitions; ++i) {
            auto it = std::find(partition[i].begin(), partition[i].end(), vertex_to_stay);
            if (it != partition[i].end()) {
                partition[i].erase(it);
                break;
            }
        }

        thread_partitions[neighbor_partition - i1].push_back(vertex_to_move);
        thread_partitions[vertex_partition - i1].push_back(vertex_to_stay);

        // Update the gains for the affected vertices
        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex_to_move) {
                if (std::find(partition[neighbor_partition].begin(), partition[neighbor_partition].end(), neighbor) != partition[neighbor_partition].end()) {
                    vertex_gains[neighbor] -= 2.0 * coarsedGraph.getEdgeWeight(vertex_to_move, neighbor);
                }
                else {
                    vertex_gains[neighbor] += 2.0 * coarsedGraph.getEdgeWeight(vertex_to_move, neighbor);
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

    // Create a list to store the partitioning for each thread.
    std::vector<std::vector<int>> thread_partitions(nthreads);

    // Copy the initial_partitions for each thread.
    std::copy(initial_partitions.begin(), initial_partitions.end(), thread_partitions.begin());

    // Create and start the threads.
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        int start_vertex = i * chunk_size;
        int end_vertex = (i + 1) * chunk_size;
        if (i == nthreads - 1) {
            end_vertex = graph.numVertices();
        }

        threads.emplace_back([i, start_vertex, end_vertex, num_iterations, &thread_partitions]() {
            RefinementStep(thread_partitions, start_vertex, end_vertex, num_iterations, thread_partitions);
        });
    }

    // Wait for all threads to finish.
    for (auto& thread : threads) {
        thread.join();
    }

    // Combine the partitioning results from all threads.
    initial_partitions.clear();
    for (const auto& partition : thread_partitions) {
        initial_partitions.emplace_back(partition);
    }
}

void WriteOutputToFile(const std::vector<std::vector<int>>& partitions, string outputfile){
    //TO DO
    int partitionIndex = 0;
    for (const auto& partition : partitions) {
        std::cout << "Partition " << partitionIndex << ": ";
        for (const auto& vertex : partition) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
        partitionIndex++;
    }
}


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile, string outputfile){
    LoadGraphFromMemory(inputfile);      //load the graph from file

    coarsedGraph = Coarsening(graph);        //Coarse the initial graph

    std::vector<std::vector<int>> initial_partitions = InitialPartitioning(npartitions);

    std::vector<std::vector<int>> uncoarsened_partitions = UncoarsePartitions(coarsedGraph, initial_partitions);

    Graph restoredGraph = Uncoarsening(coarsedGraph);

    /*
    MultithreadedRefinement(nthreads, initial_partitions);

    cout << "Refinement done" << endl;

    // Uncoarsen the partitioned graph to the original size
    std::vector<std::vector<int>> final_partitions = Uncoarsening(initial_partitions);

    cout << "Final partitions created" << endl;

    // Write the final partitioning to an output file (optional)
    WriteOutputToFile(final_partitions, outputfile);

    cout << "Output wrote in the file" << endl;
     */
}


int main() {

    // Add vertices and edges to the global graph instance
    graph.addVertex(0, 1.0);
    graph.addVertex(1, 2.0);
    graph.addVertex(2, 3.0);
    graph.addVertex(3, 4.0);
    graph.addVertex(4, 5.0);
    graph.addVertex(5, 6.0);
    graph.addVertex(6, 7.0);
    graph.addVertex(7, 8.0);
    graph.addVertex(8, 9.0);
    graph.addVertex(9, 10.0);
    graph.addVertex(10, 11.0);
    graph.addVertex(11, 12.0);
    graph.addVertex(12, 13.0);
    graph.addVertex(13, 14.0);
    graph.addVertex(14, 15.0);

    graph.addEdge(0, 1, 0.5);
    graph.addEdge(1, 2, 0.7);
    graph.addEdge(2, 0, 0.3);
    graph.addEdge(0, 3, 0.2);
    graph.addEdge(1, 3, 0.4);
    graph.addEdge(2, 3, 0.1);
    graph.addEdge(3, 4, 0.8);
    graph.addEdge(0, 4, 0.6);
    graph.addEdge(1, 4, 0.9);
    graph.addEdge(4, 5, 0.3);
    graph.addEdge(5, 6, 0.1);
    graph.addEdge(6, 7, 0.4);
    graph.addEdge(7, 8, 0.6);
    graph.addEdge(8, 9, 0.2);
    graph.addEdge(9, 10, 0.5);
    graph.addEdge(10, 11, 0.7);
    graph.addEdge(11, 12, 0.9);
    graph.addEdge(12, 13, 0.8);
    graph.addEdge(13, 14, 0.3);
    graph.addEdge(14, 0, 0.4);

    // Set the number of threads and partitions
    int nthreads = 4; // Change this to the desired number of threads
    int npartitions = 3; // Change this to the desired number of partitions
    float maxdeviation = 1.1; // Change this to the desired max deviation
    std::string inputfile = "input_graph.txt"; // Change this to the input file name
    std::string outputfile = "output_partition.txt"; // Change this to the output file name

    // Call the algorithm function
    MultithreadedMETIS(nthreads, npartitions, maxdeviation, inputfile, outputfile);

    return 0;
}
