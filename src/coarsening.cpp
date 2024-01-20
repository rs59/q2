#include "graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <condition_variable>

int coarsestGraphSize = 50;

int current_ctr;

// Function to compute matching between vertices for coarsening
// This function calculates a matching between vertices in the graph, ensuring that each vertex is matched with an unmatched neighbor.
void ComputeMatchingSubset(Graph& graph_cm, std::unordered_map<int, int>& matching, std::unordered_map<int, int>& matching_to_collapse, std::unordered_set<int>& unmatchedVertices, std::unordered_set<int>& isolatedVertices,  std::mutex& mutex, std::mutex&barrierMutex,
                           const std::unordered_map<int, double>& vertices, int start, int end, int numThreads, std::condition_variable& cv) {

    //std::cout << "Thread[" << start << "," << end << "] started" << std::endl;
    mutex.unlock();  //let an other thread start

    // Create a temporary set to store unmatched vertices
    std::unordered_set<int> tempUnmatchedVertices;

    // Populate the set of unmatched vertices without locking
    auto vertexIterator = vertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != vertices.end(); ++i, ++vertexIterator) {
        tempUnmatchedVertices.insert(vertexIterator->first);
    }

    // Lock the mutex and update the main unmatchedVertices set
    mutex.lock();
    unmatchedVertices.insert(tempUnmatchedVertices.begin(), tempUnmatchedVertices.end());
    mutex.unlock();

    //barrier implementation with cv
    std::unique_lock<std::mutex> lock(barrierMutex);
    current_ctr++;
    if (current_ctr == numThreads) {
        // All threads have arrived, notify all waiting threads
        current_ctr = 0;
        //std::cout << "Last thread arrived at the barrier, calling notify all" << std::endl;
        cv.notify_all();
    } else {
        // Not all threads have arrived, wait
        cv.wait(lock, [] { return current_ctr == 0; });
    }

    vertexIterator = vertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != vertices.end(); ++i, ++vertexIterator) {
        //Acquire the lock
        std::unique_lock<std::mutex> lock(mutex);

        int vertexId = vertexIterator->first;

        //std::cout << "Checking node " << vertexId << std::endl;

        if (matching.find(vertexId) != matching.end()) {
            //Release the lock
            //std::cout << "Already matched vertex, skipping" << std::endl;
            lock.unlock();
            continue; // Skip already matched vertices
        }

        const std::vector<int>& neighbors = graph_cm.getNeighbors(vertexId);
        int matched_neighbor = -1;

        if(neighbors.size() == 0){
            //isolated node
            isolatedVertices.insert(vertexId);
        }

        // Find the first unmatched neighbor
        for (int neighbor : neighbors) {
            //std::cout << "Neigbor: " << neighbor << std::endl;
            if (unmatchedVertices.find(neighbor) != unmatchedVertices.end()) {
                //std::cout << "Matched{" << vertexId << "," << neighbor << "}" << std::endl;
                matched_neighbor = neighbor;
                break;
            }
        }


        //if one of the neighbors is unmatched, match the evaluated vertex with that one and update data structures
        if (matched_neighbor != -1) {
            //std::cout << "Updating matching data structures" << std::endl;
            matching[vertexId] = matched_neighbor;
            matching[matched_neighbor] = vertexId;
            matching_to_collapse[vertexId] = matched_neighbor;
            unmatchedVertices.erase(vertexId);
            unmatchedVertices.erase(matched_neighbor);
        }

        //Release the lock
        lock.unlock();
    }
}

// Function to compute matching between vertices for coarsening
//This function divide the vertices among the different nthreads and call the function that matches this subset
std::unordered_map<int, int> ComputeMatching(Graph& graph_cm, int numThreads, std::unordered_map<int, int>& matching_to_collapse, std::unordered_set<int>& unmatchedVertices) {
    std::unordered_map<int, int> matching;

    std::unordered_set<int> isolatedVertices;
    std::unordered_map<int, double> vertices = graph_cm.getVertices();

    std::cout << "ComputeMatching entered" << std::endl;

    // Create a single mutex for synchronization
    std::mutex mutex;

    //Create a barrier mutex
    std::mutex barrierMutex;

    // Divide the work among three threads
    int verticesPerThread = vertices.size() / numThreads;

    current_ctr = 0;  //reset barrier ctr
    std::condition_variable cv;

    // Create and join threads in a loop
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; i++) {
        mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
        int start = i * verticesPerThread;
        int end = (i + 1) * verticesPerThread;

        if (i == numThreads - 1) {
            end = vertices.size(); // Handle the last thread's range
        }

        threads.emplace_back([&] {
            ComputeMatchingSubset(graph_cm, matching, matching_to_collapse, unmatchedVertices, isolatedVertices, std::ref(mutex), std::ref(barrierMutex), vertices, start, end, numThreads, std::ref(cv));
        });
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }

    //Match isolated nodes
    auto it = isolatedVertices.begin();

    // Iterate until we reach the end of isolatedVertices or just before the last element
    while (it != isolatedVertices.end() && std::next(it) != isolatedVertices.end()) {
        // Match the current node with the next node
        int node1 = *it;
        int node2 = *std::next(it);

        // Store the matching pair in the map
        matching[node1] = node2;
        matching[node2] = node1;
        matching_to_collapse[node1] = node2;
        unmatchedVertices.erase(node1);
        unmatchedVertices.erase(node2);

        // Move the iterator to the next pair of isolated vertices
        std::advance(it, 2);
    }

    std::cout << "ComputeMatching exited" << std::endl;

    return matching;
}

// Function to collapse vertices based on matching information obtained in ComputeMatching
// This function collapses vertices in the graph based on the matching information obtained during coarsening, adding the collapsed vertices (that are less than original ones)
// to the coarsed graph
//matching = matching found in ComputeMatching function, used to collapse nodes
//vertices_map = originally empty, an entry with index equal to the vertex will be created each time a vertex is collapsed and added to coarsed_graph, containing the value of the new node's value
//Example: node 1 is matcehd with node 3 and will be collapsed in new node 5, so vertices_map[1] == vertices_map[3] == 5;
void CollapseVerticesSubset(Graph& graph_cv, Graph& coarsed_graph, std::unordered_map<int, int>& matching, std::unordered_map<int, int>& matching_to_collapse, std::unordered_set<int>& unmatchedVertices, std::mutex& mutex, std::unordered_map<int,int>& vertices_map, std::unordered_map<int, double>& vertices, int& newVertexCtr, int start, int end) {
    //Release the lock so that other threads can be created
    //std::cout << "Thread[" << start << "," << end << "] started" << std::endl;
    mutex.unlock();

    auto vertexIterator = vertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != vertices.end(); ++i, ++vertexIterator) {

        int vertexId = vertexIterator->first;

        //std::cout << "Evaluating vertex: " << vertexId << std::endl;

        auto matched_vertex_it = matching_to_collapse.find(vertexId);

        if(matched_vertex_it != matching_to_collapse.end() || unmatchedVertices.find(vertexId) != unmatchedVertices.end()){
            if (matched_vertex_it != matching_to_collapse.end()) {
                // Key (vertexID) was found in the map
                int matched_neighbor = matched_vertex_it->second; // Get the matched neighbor
                double new_weight = graph_cv.getVertexWeight(vertexId) + graph_cv.getVertexWeight(matched_neighbor);  //Sum the weights of the collapsed vertices
                mutex.lock();
                coarsed_graph.addVertex(newVertexCtr, new_weight);  //collapse the vertices

                //Associate to both the vertices in vertices_map the new value of the node, and increase it for the next couple
                vertices_map[vertexId] = newVertexCtr;
                vertices_map[matched_neighbor] = newVertexCtr;
                newVertexCtr = newVertexCtr + 1;
                mutex.unlock();
            } else {
                mutex.lock();
                //unmatched vertex
                coarsed_graph.addVertex(newVertexCtr, graph_cv.getVertexWeight(vertexId));
                vertices_map[vertexId] = newVertexCtr;
                newVertexCtr = newVertexCtr + 1;
                mutex.unlock();
            }
        }
    }
}

// Function to collapse vertices based on matching information
// This function collapses vertices in the graph based on the matching information obtained during coarsening.
void CollapseVertices(Graph& coarsed_graph, Graph& graph_cv, std::unordered_map<int, int>& matching, std::unordered_map<int, int>& matching_to_collapse, std::unordered_set<int>& unmatchedVertices, int numThreads) {

    std::cout << "CollapseVertices entered" << std::endl;

/*
    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();
    //save in the new coarsed_graph the previous coarsening mappings
    for(int i = 0; i < graph_cv.getCoarsingLevel(); i++){
        coarsed_graph.pushBackMapping(graph_cv.getMapping(i));
        coarsed_graph.pushBackMappingEdges(graph_cv.getMappingEdges(i));
        coarsed_graph.pushBackVerticesWeights(graph_cv.getMappingVerticesWeights(i));
    }
    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Restore mappings: " << seconds << " seconds" << std::endl;
*/

    std::unordered_map<int,int> vertices_map;  //store the mapping between coarsed vertices and original ones {original vertex, new vertex}
    std::unordered_map<int, double> vertices = graph_cv.getVertices();

    // Create a single mutex for synchronization
    std::mutex mutex;

    // Divide the work among three threads
    int verticesPerThread = vertices.size() / numThreads;

    int newVertexCtr = 0;  //this will be shared among threads to keep count of new vertices

    // Create and join threads in a loop
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; i++) {
        mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
        int start = i * verticesPerThread;
        int end = (i + 1) * verticesPerThread;

        if (i == numThreads - 1) {
            end = vertices.size(); // Handle the last thread's range
        }

        threads.emplace_back([&] {
            CollapseVerticesSubset(graph_cv, coarsed_graph, matching, matching_to_collapse, unmatchedVertices, std::ref(mutex), vertices_map, vertices, newVertexCtr, start, end);
        });
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }

    coarsed_graph.pushBackMapping(vertices_map);
    coarsed_graph.pushBackMappingEdges(graph_cv.getEdgeWeights());
    coarsed_graph.pushBackVerticesWeights(graph_cv.getVertices());

    std::cout << "CollapseVertices exited" << std::endl;
}

// Function to update the edge weights after collapsing vertices for a subset of edges
// This function is executed by each thread to process a subset of edges.
void UpdateEdgeWeightsSubset(Graph& graph_ue, const std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights, std::unordered_map<int,int>& vertices_map, std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights_updated, std::mutex& mutex, int start, int end) {
    //std::cout << "Thread[" << start << "," << end << "]" << std::endl;
    mutex.unlock();

    std::unordered_map<std::pair<int, int>, double, HashPair> tempEdge_weights;

    auto edgeIterator = edge_weights.begin();
    std::advance(edgeIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && edgeIterator != edge_weights.end(); ++i, ++edgeIterator) {
        int nodeA = edgeIterator->first.first;
        int nodeB = edgeIterator->first.second;
        double value = edgeIterator->second;

        if (vertices_map.find(nodeA) != vertices_map.end() && vertices_map.find(nodeB) != vertices_map.end()) {
            //Find the new value of the 2 nodes after having collapsed them, using vertices_map data structure
            int coarsedSource = vertices_map.at(nodeA);
            int coarsedDestination = vertices_map.at(nodeB);

            // Skip processing if coarsedSource and coarsedDestination are the same
            if (coarsedSource == coarsedDestination) {
                continue;
            }

            tempEdge_weights.emplace(std::make_pair(coarsedSource,coarsedDestination), value);
        }
    }

    int coarsedSource;
    int coarsedDestination;
    double value;
    mutex.lock();
    for (auto it = tempEdge_weights.begin(); it != tempEdge_weights.end(); ++it) {
        coarsedSource = it->first.first;
        coarsedDestination = it->first.second;
        value = it->second;

        // Check if an edge already exists between the coarsed vertices
        if (graph_ue.containsEdge(coarsedSource, coarsedDestination)) {
            double existingWeight = graph_ue.getEdgeWeight(coarsedSource, coarsedDestination);
            graph_ue.addEdge(coarsedSource, coarsedDestination, existingWeight + (value/2)); //value is /2 because we have double edges
        } else {
            graph_ue.addEdge(coarsedSource, coarsedDestination, (value/2)); //value is /2 because there are double edges
        }
    }
    mutex.unlock();
}

// Function to update edge weights after collapsing vertices
// This function updates edge weights in the graph after collapsing vertices.
//After this function, the collapsed graph will have both collapsed nodes and updated edges between them
void UpdateEdgeWeights(Graph& graph_ue, const std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights, int numThreads) {
    std::cout << "UpdateEdgeWeights entered" << std::endl;

    std::unordered_map<int,int> vertices_map = graph_ue.getMapping(graph_ue.getCoarsingLevel() - 1);

    std::unordered_map<std::pair<int, int>, double, HashPair> edge_weights_updated;

    // Create a single mutex for synchronization
    std::mutex mutex;

    // Divide the work among three threads
    int verticesPerThread = edge_weights.size() / numThreads;

    // Create and join threads in a loop
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; i++) {
        mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
        int start = i * verticesPerThread;
        int end = (i + 1) * verticesPerThread;

        if (i == numThreads - 1) {
            end = edge_weights.size(); // Handle the last thread's range
        }

        threads.emplace_back([&] {
            UpdateEdgeWeightsSubset(graph_ue, edge_weights, vertices_map, edge_weights_updated, std::ref(mutex), start, end);
        });
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }
    std::cout << "UpdateEdgeWeights exited" << std::endl;
}

Graph Coarsening(Graph& graph, int nthreads, int npartitions){
    Graph coarsened_graph = graph;

    //Calc the number of isolated nodes
    int isolatedNodes = 0;

    for (const auto& entry : graph.getVertices()) {
        int vertexId = entry.first;    // Get the key (int)

        if(graph.getNeighbors(vertexId).size() == 0){
            isolatedNodes++;
        }
    }

    int coarsest_graph_size = 0;
    if(((graph.size()/3) + isolatedNodes) < coarsestGraphSize && ((graph.size()/3) + isolatedNodes) >= npartitions){
        coarsest_graph_size = (graph.size()/3) + isolatedNodes;
    }else{
        coarsest_graph_size = coarsestGraphSize;
    }

    int i = 0;

    std::cout << "Starting coarsening phase" << std::endl;

    while(coarsened_graph.size() > coarsest_graph_size){
        std::cout << "New iteration" << std::endl;
        Graph temp_graph = coarsened_graph;
        std::cout << "Graph size: " << temp_graph.getVertices().size() << std::endl;
        // Start the clock
        auto start_time = std::chrono::high_resolution_clock::now();
        std::unordered_map<int, int> matching_to_collapse;
        std::unordered_set<int> unmatchedVertices;
        std::unordered_map<int, int> matching = ComputeMatching(temp_graph, nthreads, matching_to_collapse, unmatchedVertices);
        // Stop the clock
        auto end_time = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        // Convert the duration to a double value in seconds
        double seconds = duration.count() / 1e6;
        // Print the execution time
        std::cout << "Compute matching: " << seconds << " seconds" << std::endl;

        std::cout << "Matching to collapse size: " << matching_to_collapse.size() << ", unmatched nodes: " << unmatchedVertices.size() << ", TOT: " << ((matching_to_collapse.size()*2) + unmatchedVertices.size()) << std::endl;

        // Start the clock
        start_time = std::chrono::high_resolution_clock::now();
        coarsened_graph.clearButMappings();
        CollapseVertices(coarsened_graph, temp_graph, matching, matching_to_collapse, unmatchedVertices, nthreads);   //graph with collapsed vertices and no edges
        // Stop the clock
        end_time = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        // Convert the duration to a double value in seconds
        seconds = duration.count() / 1e6;
        // Print the execution time
        std::cout << "Collapse vertices: " << seconds << " seconds" << std::endl;

        // Start the clock
        start_time = std::chrono::high_resolution_clock::now();
        std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights = temp_graph.getEdgeWeights();
        UpdateEdgeWeights(coarsened_graph, edgeWeights, nthreads);
        // Stop the clock
        end_time = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        // Convert the duration to a double value in seconds
        seconds = duration.count() / 1e6;
        // Print the execution time
        std::cout << "Update edge weights: " << seconds << " seconds" << std::endl;

        i++;
        std::cout << "We reached the end of the cycle, new size: " << coarsened_graph.size() << ", condition size: " << coarsest_graph_size << std::endl;
    }

    return coarsened_graph;

}