#include<iostream>
#include<graph.h>
#include<vector>
#include<algorithm>

int NUM_ITERATIONS = 100; //TO BE DEFINED
int coarsest_graph_size = 10 //TO BE DEFINED

//To be defined in LoadGraphFromMemory 
Graph graph;
Graph coarsedGraph;

void LoadGraphFromMemory(string inputfile){
    //TO-DO
}

// Function to compute matching between vertices for coarsening
std::vector<std::pair<int, int>> ComputeMatching(const Graph& graph_cm) {
    std::vector<std::pair<int, int>> matching;
    std::vector<bool> is_matched(graph_cm.size(), false);

    for (int vertex = 0; vertex < graph_cm.size(); vertex++) {
        if (is_matched[vertex]) {
            continue; // Skip already matched vertices
        }

        const std::vector<int>& neighbors = graph_cm[vertex];
        int matched_neighbor = -1;

        // Find the first unmatched neighbor
        for (int neighbor : neighbors) {
            if (!is_matched[neighbor]) {
                matched_neighbor = neighbor;
                break;
            }
        }

        if (matched_neighbor != -1) {
            matching.emplace_back(vertex, matched_neighbor);
            is_matched[vertex] = true;
            is_matched[matched_neighbor] = true;
        }
    }

    return matching;
}

// Function to collapse matched vertices into a coarser representation
Graph CollapseVertices(const Graph& graph_cv, const std::vector<std::pair<int, int>>& matching) {
    Graph coarse_graph;
    std::unordered_set<int> merged_vertices;

    for (int vertex = 0; vertex < graph_cv.size(); vertex++) {
        if (merged_vertices.find(vertex) == merged_vertices.end()) {
            auto matched_vertex_it = std::find_if(matching.begin(), matching.end(), [vertex](const std::pair<int, int>& match) {
                return match.first == vertex || match.second == vertex;
            });

            if (matched_vertex_it == matching.end()) {
                coarse_graph.push_back(graph_cv[vertex]);
            } else {
                int matched_vertex = (vertex == matched_vertex_it->first) ? matched_vertex_it->second : matched_vertex_it->first;
                merged_vertices.insert(matched_vertex);
                merged_vertices.insert(vertex);
                std::vector<int> merged_neighbors = graph_cv[vertex];
                merged_neighbors.insert(merged_neighbors.end(), graph_cv[matched_vertex].begin(), graph_cv[matched_vertex].end());
                coarse_graph.push_back(merged_neighbors);
            }
        }
    }

    return coarse_graph;
}

// Function to update the edge weights after collapsing vertices
std::vector<double> UpdateEdgeWeights(const Graph& graph_ue, const std::vector<std::pair<int, int>>& matching, const std::vector<double>& edge_weights) {
    std::vector<double> coarse_edge_weights(graph_ue.size(), 0.0);

    for (int vertex = 0; vertex < graph_ue.size(); vertex++) {
        int matched_vertex = -1;
        for (const auto& match : matching) {
            if (match.first == vertex || match.second == vertex) {
                matched_vertex = (vertex == match.first) ? match.second : match.first;
                break;
            }
        }

        if (matched_vertex == -1) {
            coarse_edge_weights[vertex] = edge_weights[vertex];
        } else {
            coarse_edge_weights[vertex] = edge_weights[vertex] + edge_weights[matched_vertex];
            // Reset the edge weight of the matched vertex since it's merged with another vertex
            edge_weights[matched_vertex] = 0.0;
        }
    }

    return coarse_edge_weights;
}

// Function to update the vertex mapping after collapsing vertices
std::vector<int> UpdateVertexMapping(const std::vector<std::pair<int, int>>& matching, const std::vector<int>& vertex_mapping) {
    std::vector<int> new_vertex_mapping(vertex_mapping.size(), -1);

    for (int vertex = 0; vertex < vertex_mapping.size(); vertex++) {
        int matched_vertex = -1;
        for (const auto& match : matching) {
            if (match.first == vertex) {
                matched_vertex = match.second;
                break;
            } else if (match.second == vertex) {
                matched_vertex = match.first;
                break;
            }
        }

        if (matched_vertex != -1) {
            int coarser_vertex = vertex_mapping[matched_vertex];
            if (new_vertex_mapping[coarser_vertex] != -1) {
                new_vertex_mapping[vertex] = new_vertex_mapping[coarser_vertex];
            } else {
                int new_index = new_vertex_mapping.size(); // Next available index
                new_vertex_mapping[coarser_vertex] = new_index;
                new_vertex_mapping[vertex] = new_index;
            }
        } else {
            // If the current vertex is not matched, keep its existing mapping
            new_vertex_mapping[vertex] = vertex_mapping[vertex];
        }
    }

    return new_vertex_mapping;
}

Graph Coarsening() {
    Graph coarse_graph = initial_graph; // Copy the initial graph representation
    std::vector<std::vector<double>> edge_weights = graph.getEdgeWeights();

    // Initialize a vector to store the mapping of vertices from finer to coarser levels
    std::vector<int> vertex_mapping(initial_graph.size());
    for (int i = 0; i < initial_graph.size(); i++) {
        vertex_mapping[i] = i; // Initialize the mapping to be identity at the beginning
    }

    // Start the coarsening process until the coarsest level is reached
    while (coarse_graph.size() > coarsest_graph_size) {
        // Compute a matching of vertices in the current level
        std::vector<std::pair<int, int>> matching = ComputeMatching(coarse_graph); // Compute a matching of vertices in the coarse_graph

        // Collapse the matched vertices to form the coarser graph
        coarse_graph = CollapseVertices(coarse_graph, matching); // Combine matched vertices to form a coarser graph

        // Update the edge weights based on the collapsed vertices
        UpdateEdgeWeights(coarse_graph, matching, edge_weights); // Update edge weights based on the matching

        // Update the vertex mapping for the finer-to-coarser level mapping
        vertex_mapping = UpdateVertexMapping(matching, vertex_mapping); // Update the vertex mapping for the next level

        // Repeat the coarsening process until the coarsest level is reached
    }

    // Return the coarsest level graph
    return coarse_graph;
}

std::vector<std::vector<int>> InitialPartitioning(int npartitions) {
    // Calculate the total weight of all vertices in the graph
    double total_weight = 0;
    for (const auto& vertexPair : graph.getVertices()) {
        total_weight += vertexPair.second;
    }

    // Calculate the target weight for each partition
    double target_weight = total_weight / npartitions;

    // Initialize the partitions and their current weight
    std::vector<std::vector<int>> partitions(npartitions);
    std::vector<double> partition_weights(npartitions, 0);

    // Create an index vector to sort vertices based on their weights
    std::vector<int> sorted_vertices;
    for (const auto& vertexPair : graph.getVertices()) {
        sorted_vertices.push_back(vertexPair.first);
    }

    // Sort vertices in descending order of their weights
    std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&graph](int a, int b) {
        return graph.getVertexWeight(a) > graph.getVertexWeight(b);
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
        partition_weights[min_weight_partition] += graph.getVertexWeight(vertex);
    }

    return partitions;
}

// Function to perform the refinement process in a multithreaded manner.
void MultithreadedRefinement(int nthreads, std::vector<std::vector<int>>& initial_partitions) {
    int num_iterations = NUM_ITERATIONS;

    // Calculate the number of vertices in each chunk for each thread.
    int chunk_size = graph.getVertices().size() / nthreads;

    // Create a list to store the partitioning for each thread.
    std::vector<std::vector<int>> thread_partitions(nthreads);

    // Copy the initial_partitions for each thread.
    std::copy(initial_partitions.begin(), initial_partitions.end(), thread_partitions.begin());

    // Create the mutex to synchronize the modification of thread_partitions.
    std::mutex mutex;

    // Create and start the threads.
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        int start_vertex = i * chunk_size;
        int end_vertex = (i + 1) * chunk_size;
        if (i == nthreads - 1) {
            end_vertex = graph.getVertices().size();
        }
        threads.emplace_back(RefinementStep, i, start_vertex, end_vertex, num_iterations, std::ref(thread_partitions));
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

// Function to perform a single step of the refinement process using the Kernighan-Lin algorithm
// to improve the partition quality.
void RefinementStep(std::vector<std::vector<int>>& partition) {
    
    std::vector<std::vector<double>> edge_weights = graph.getEdgeWeights();

    // Get the number of vertices and partitions
    int num_vertices = graph.getVertices().size();
    int num_partitions = partition.size();

    // Calculate the target weight for each partition
    double total_weight = 0.0;
    for (const auto& entry : edgeWeights) {
        total_weight += entry.second;
    }
    double target_weight = total_weight / num_partitions;

    // Initialize the partition gain for each vertex
    std::vector<double> vertex_gains(num_vertices, 0.0);

    // Calculate the initial gain for each vertex in its current partition
    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        int vertex_partition = -1;
        for (int i = 0; i < num_partitions; ++i) {
            if (std::find(partition[i].begin(), partition[i].end(), vertex) != partition[i].end()) {
                vertex_partition = i;
                break;
            }
        }

        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex) {
                auto edge_weight_it = edgeWeights.find({vertex, neighbor});
                double weight = (edge_weight_it != edgeWeights.end()) ? edge_weight_it->second : 0.0;

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
    for (int iter = 0; iter < num_vertices / 2; ++iter) { // A heuristic to limit the number of swaps

        // Find the best pair of vertices to swap
        double best_gain = 0.0;
        int vertex_to_move = -1;
        int vertex_to_stay = -1;

        for (int vertex = 0; vertex < num_vertices; ++vertex) {
            if (locked_vertices.find(vertex) == locked_vertices.end() && (vertex_to_move == -1 || vertex_gains[vertex] > best_gain)) {
                vertex_to_move = vertex;
                best_gain = vertex_gains[vertex];
            }
        }

        locked_vertices.insert(vertex_to_move);

        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex_to_move) {
                int neighbor_partition = -1;
                for (int i = 0; i < num_partitions; ++i) {
                    if (std::find(partition[i].begin(), partition[i].end(), neighbor) != partition[i].end()) {
                        neighbor_partition = i;
                        break;
                    }
                }

                double gain = vertex_gains[vertex_to_move];
                if (neighbor_partition != -1) {
                    gain -= 2.0 * edge_weights[vertex_to_move][neighbor];
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

        partition[neighbor_partition].push_back(vertex_to_move);
        partition[vertex_partition].push_back(vertex_to_stay);

        // Update the gains for the affected vertices
        for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
            if (neighbor != vertex_to_move) {
                if (std::find(partition[neighbor_partition].begin(), partition[neighbor_partition].end(), neighbor) != partition[neighbor_partition].end()) {
                    vertex_gains[neighbor] -= 2.0 * edge_weights[vertex_to_move][neighbor];
                }
                else {
                    vertex_gains[neighbor] += 2.0 * edge_weights[vertex_to_move][neighbor];
                }
            }
        }
    }
}


void MultithreadedMETIS(int nthreads, int npartitions, float maxdeviation, string inputfile){
    LoadGraphFromMemory(inputfile);      //load the graph from file
    coarsedGraph = Coarsening();        //Coarse the initial graph

    std::vector<std::vector<int>> initial_partitions = InitialPartitioning(npartitions);

    MultithreadedRefinement(nthreads, initial_partitions);

    WriteOutputToFile(final_partitions, outputfile);
}