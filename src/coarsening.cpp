#include "graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <condition_variable>

int coarsestGraphSize = 30000;

int current_ctr;

// Function to compute matching between vertices for coarsening
void ComputeMatchingSubset(Graph& graph_cm, std::unordered_map<int, int>& matching, std::unordered_set<int>& unmatchedVertices,  std::mutex& mutex, std::mutex&barrierMutex,
                           const std::unordered_map<int, double>& vertices, int start, int end, int numThreads, std::condition_variable& cv) {

    //std::cout << "Thread[" << start << "," << end << "] started" << std::endl;
    mutex.unlock();  //let an other thread start

    //Populate the set of unmatched vertices
    auto vertexIterator = vertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != vertices.end(); ++i, ++vertexIterator) {
        mutex.lock();
        unmatchedVertices.insert(vertexIterator->first);
        mutex.unlock();
    }

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

        // Find the first unmatched neighbor
        for (int neighbor : neighbors) {
            //std::cout << "Neigbor: " << neighbor << std::endl;
            if (unmatchedVertices.find(neighbor) != unmatchedVertices.end()) {
                //std::cout << "Matched{" << vertexId << "," << neighbor << "}" << std::endl;
                matched_neighbor = neighbor;
                break;
            }
        }



        if (matched_neighbor != -1) {
            //std::cout << "Updating matching data structures" << std::endl;
            matching[vertexId] = matched_neighbor;
            matching[matched_neighbor] = vertexId;
            unmatchedVertices.erase(vertexId);
            unmatchedVertices.erase(matched_neighbor);
        }

        //Release the lock
        lock.unlock();
    }
}

// Function to compute matching between vertices for coarsening
std::unordered_map<int, int> ComputeMatching(Graph& graph_cm, int numThreads) {
    std::unordered_map<int, int> matching;
    std::unordered_set<int> unmatchedVertices;
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
            ComputeMatchingSubset(graph_cm, matching, unmatchedVertices, std::ref(mutex), std::ref(barrierMutex), vertices, start, end, numThreads, std::ref(cv));
        });
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }
    std::cout << "ComputeMatching exited" << std::endl;

    return matching;
}


void CollapseVerticesSubset(Graph& graph_cv, Graph& coarsed_graph, std::unordered_map<int, int>& matching, std::mutex& mutex, std::unordered_map<int,int>& vertices_map, std::unordered_map<int, double>& vertices, int& newVertexCtr, int start, int end) {
    //Release the lock so that other threads can be created
    //std::cout << "Thread[" << start << "," << end << "] started" << std::endl;
    mutex.unlock();

    auto vertexIterator = vertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != vertices.end(); ++i, ++vertexIterator) {
        //Acquire the lock
        std::unique_lock<std::mutex> lock(mutex);

        int vertexId = vertexIterator->first;

        //std::cout << "Evaluating vertex: " << vertexId << std::endl;

        if (vertices_map.find(vertexId) != vertices_map.end()) {
            //Release the lock
            //std::cout << "Already inserted vertex, skipping" << std::endl;
            lock.unlock();
            continue; // Skip already inserted vertices
        }
        auto matched_vertex_it = matching.find(vertexId);

        if (matched_vertex_it != matching.end()) {
            // Key (vertexID) was found in the map
            //std::cout << "Collapsing vertices " << vertexId << " and " << matched_vertex_it->second << std::endl;
            int matched_neighbor = matched_vertex_it->second; // Get the matched neighbor
            double new_weight = graph_cv.getVertexWeight(vertexId) + graph_cv.getVertexWeight(matched_neighbor);
            coarsed_graph.addVertex(newVertexCtr, new_weight);  //collapse the vertices
            vertices_map[vertexId] = newVertexCtr;
            vertices_map[matched_neighbor] = newVertexCtr;
            newVertexCtr = newVertexCtr + 1;
        } else {
            // Key (vertexID) was not found in the map
            //std::cout << "Adding vertex " << vertexId << " without collapsing" << std::endl;
            coarsed_graph.addVertex(newVertexCtr, graph_cv.getVertexWeight(vertexId));
            vertices_map[vertexId] = newVertexCtr;
            newVertexCtr = newVertexCtr + 1;
        }
        //Release the lock
        lock.unlock();
    }
}

Graph CollapseVertices(Graph& graph_cv, std::unordered_map<int, int>& matching, int numThreads) {
    Graph coarsed_graph;

    std::cout << "CollapseVertices entered" << std::endl;

    //save in the new coarsed_graph the previous coarsening mappings
    for(int i = 0; i < graph_cv.getCoarsingLevel(); i++){
        coarsed_graph.pushBackMapping(graph_cv.getMapping(i));
        coarsed_graph.pushBackMappingEdges(graph_cv.getMappingEdges(i));
        coarsed_graph.pushBackVerticesWeights(graph_cv.getMappingVerticesWeights(i));
    }

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
            CollapseVerticesSubset(graph_cv, coarsed_graph, matching, std::ref(mutex), vertices_map, vertices, newVertexCtr, start, end);
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

    return coarsed_graph;
}

// Function to update the edge weights after collapsing vertices
void UpdateEdgeWeightsSubset(Graph& graph_ue, const std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights, std::unordered_map<int,int>& vertices_map, std::unordered_map<std::pair<int, int>, double, HashPair>& edge_weights_updated, std::mutex& mutex, int start, int end) {
    //std::cout << "Thread[" << start << "," << end << "]" << std::endl;
    mutex.unlock();

    auto edgeIterator = edge_weights.begin();
    std::advance(edgeIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && edgeIterator != edge_weights.end(); ++i, ++edgeIterator) {
        int nodeA = edgeIterator->first.first;
        int nodeB = edgeIterator->first.second;
        double value = edgeIterator->second;

        if (vertices_map.find(nodeA) != vertices_map.end() && vertices_map.find(nodeB) != vertices_map.end()) {
            int coarsedSource = vertices_map.at(nodeA);
            int coarsedDestination = vertices_map.at(nodeB);

            // Skip processing if coarsedSource and coarsedDestination are the same
            if (coarsedSource == coarsedDestination) {
                continue;
            }

            // Check if an edge already exists between the coarsed vertices
            if (graph_ue.containsEdge(coarsedSource, coarsedDestination)) {
                double existingWeight = graph_ue.getEdgeWeight(coarsedSource, coarsedDestination);
                // Sum the weights if an edge already exists
                mutex.lock();
                graph_ue.addEdge(coarsedSource, coarsedDestination, existingWeight + (value/2)); //value is /2 because we have double edges
                mutex.unlock();
            } else {
                mutex.lock();
                graph_ue.addEdge(coarsedSource, coarsedDestination, (value/2));
                mutex.unlock();
            }
        }
    }
}

// Function to update the edge weights after collapsing vertices
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

Graph Coarsening(Graph& graph, int nthreads){
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
    if(((graph.size()/3) + isolatedNodes) < coarsestGraphSize){
        coarsest_graph_size = (graph.size()/3) + isolatedNodes;
    }else{
        coarsest_graph_size = coarsestGraphSize;
    }

    int i = 0;

    std::cout << "Starting coarsening phase" << std::endl;

    while(coarsened_graph.size() > coarsest_graph_size){
        std::cout << "New iteration" << std::endl;
        Graph temp_graph = coarsened_graph;
        std::unordered_map<int, int> matching = ComputeMatching(temp_graph, nthreads);
        coarsened_graph = CollapseVertices(temp_graph, matching, nthreads);   //graph with collapsed vertices and no edges
        std::unordered_map<std::pair<int, int>, double, HashPair> edgeWeights = temp_graph.getEdgeWeights();

        UpdateEdgeWeights(coarsened_graph, edgeWeights, nthreads);

        i++;
        std::cout << "We reached the end of the cycle, new size: " << coarsened_graph.size() << ", condition size: " << coarsest_graph_size << std::endl;
    }

    return coarsened_graph;

}