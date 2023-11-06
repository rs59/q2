#include <iostream>
#include "graph.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <condition_variable>

int crrnt_ctr;
int MAX_SWAPS_PER_ITERATION = 1000;

void UncoarseningGraph(Graph& uncoarsened_temp, std::unordered_map<int, int>& coarser_to_finer_mapping, std::unordered_map<std::pair<int, int>, double, HashPair>& edgesMapping, std::unordered_map<int, double>& verticesWeights, std::mutex& mutex, int start_mapping, int end_mapping, int start_edges, int end_edges, int numThreads, std::condition_variable& cv, std::mutex&barrierMutex){
    mutex.unlock();

    auto mappingIterator = coarser_to_finer_mapping.begin();
    std::advance(mappingIterator, start_mapping); // Move the iterator to the starting position
    for (int i = start_mapping; i < end_mapping && mappingIterator != coarser_to_finer_mapping.end(); i++, mappingIterator++) {
        //Acquire the lock
        std::unique_lock<std::mutex> lock(mutex);

        int original = mappingIterator -> first;

        auto vertexIt = verticesWeights.find(original);
        if(vertexIt != verticesWeights.end()){
            double weight = vertexIt->second;
            uncoarsened_temp.addVertex(original, weight);
        }
        lock.unlock();
    }

    //barrier implementation with cv
    std::unique_lock<std::mutex> lock(barrierMutex);
    crrnt_ctr++;
    if (crrnt_ctr == numThreads) {
        // All threads have arrived, notify all waiting threads
        crrnt_ctr = 0;
        cv.notify_all();
    } else {
        // Not all threads have arrived, wait
        cv.wait(lock, [] { return crrnt_ctr == 0; });
    }

    auto edgesIterator = edgesMapping.begin();
    std::advance(edgesIterator, start_edges); // Move the iterator to the starting position
    for (int i = start_edges; i < end_edges && edgesIterator != edgesMapping.end(); i++, edgesIterator++) {
        //Acquire the lock
        std::unique_lock<std::mutex> lock(mutex);

        int nodeA = edgesIterator->first.first;
        int nodeB = edgesIterator->first.second;
        double value = edgesIterator->second;

        uncoarsened_temp.addEdge(nodeA, nodeB, value);
        lock.unlock();
    }
}

void UncoarsePartitions(std::unordered_map<int, int>& coarser_to_finer_mapping, std::vector<std::vector<int>>& partitions, std::vector<std::vector<int>>& uncoarsed_partitions, std::mutex& mutex, int start, int end) {
    //std::cout << "Thread here" << std::endl;
    mutex.unlock();

    auto mappingIterator = coarser_to_finer_mapping.begin();
    std::advance(mappingIterator, start); // Move the iterator to the starting position
    int originalVertex;
    int coarserVertex;
    for (int i = start; i < end && mappingIterator != coarser_to_finer_mapping.end(); i++, mappingIterator++) {
        originalVertex = mappingIterator -> first;
        coarserVertex = mappingIterator -> second;

        //std::cout << "[" << originalVertex << "," << coarserVertex << "]" << std::endl;

        for(int i = 0; i < partitions.size(); i++){
            for(int j = 0; j < partitions[i].size(); j++){
                if(partitions[i][j] == coarserVertex){
                    mutex.lock();
                    uncoarsed_partitions[i].emplace_back(originalVertex);
                    mutex.unlock();
                }
            }
        }
    }
}

void uncoarseBoundaries(Graph& graph, std::unordered_map<int, int>& coarser_to_finer_mapping, std::vector<int>& boundaryVertices, std::vector<int>& uncoarsedBoundaryVertices, std::vector<std::vector<int>>& partitions,  std::mutex& mutex, int start, int end){
    //std::cout << "Thread[" << start << "," << end << "] started" << std::endl;
    mutex.unlock();

    auto mappingIterator = coarser_to_finer_mapping.begin();
    std::advance(mappingIterator, start); // Move the iterator to the starting position
    int originalVertex;
    int coarserVertex;
    for (int i = start; i < end && mappingIterator != coarser_to_finer_mapping.end(); i++, mappingIterator++) {
        originalVertex = mappingIterator -> first;
        coarserVertex = mappingIterator -> second;

        if(std::find(boundaryVertices.begin(), boundaryVertices.end(), coarserVertex) == boundaryVertices.end()){
            //The coarsed vertex is not boundary, so none of the original ones will be
            continue;
        }

        std::vector<int> neigbhors = graph.getNeighbors(originalVertex);

        // Get the partition of the current vertex
        int partitionVertex = -1; // Initialize to an invalid partition

        // Find the partition that contains the current vertex
        for (int i = 0; i < partitions.size(); ++i) {
            if (std::find(partitions[i].begin(), partitions[i].end(), originalVertex) != partitions[i].end()) {
                partitionVertex = i;
                break;
            }
        }

        // Initialize a flag to check if the vertex is a boundary vertex
        bool isBoundary = false;

        // Iterate through the neighbors of the current vertex
        for (int neighbor : neigbhors) {
            // Find the partition that contains the neighbor
            int partitionNeighbor = -1; // Initialize to an invalid partition
            for (int i = 0; i < partitions.size(); ++i) {
                if (std::find(partitions[i].begin(), partitions[i].end(), neighbor) != partitions[i].end()) {
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
            mutex.lock();
            uncoarsedBoundaryVertices.push_back(originalVertex);
            mutex.unlock();
        }
    }
}

void refinementStep(Graph& graph, std::vector<std::pair<int, double>>& vertexGains, std::vector<std::vector<int>>& partitions, std::vector<int>& boundaryVertices, std::mutex& mutex, int start, int end){
    mutex.unlock();

    // Calculate initial vertex gains
    // Iterate through each vertex from start to end
    auto vertexIterator = boundaryVertices.begin();
    std::advance(vertexIterator, start); // Move the iterator to the starting position
    for (int i = start; i < end && vertexIterator != boundaryVertices.end(); i++, vertexIterator++) {
        int vertexID = *vertexIterator;

        const std::vector<int>& neighbors = graph.getNeighbors(vertexID);

        // Get the partition of the current vertex
        int partitionVertex = -1; // Initialize to an invalid partition
        for (int i = 0; i < partitions.size(); ++i) {
            if (std::find(partitions[i].begin(), partitions[i].end(), vertexID) != partitions[i].end()) {
                partitionVertex = i;
                break;
            }
        }

        // Calculate the initial gain for each vertex in its current partition
        double initialGain = 0.0;
        for (int neighbor : neighbors) { // Assuming graph is an adjacency list
            double weight = graph.getEdgeWeight(vertexID, neighbor); // Assuming you have this data structure

            // Find the partition that contains the neighbor
            int partitionNeighbor = -1; // Initialize to an invalid partition
            for (int j = 0; j < partitions.size(); ++j) {
                if (std::find(partitions[j].begin(), partitions[j].end(), neighbor) != partitions[j].end()) {
                    partitionNeighbor = j;
                    break;
                }
            }

            if (partitionVertex == partitionNeighbor) {
                initialGain -= weight;
            } else {
                initialGain += weight;
            }
        }
        mutex.lock();
        vertexGains.emplace_back(vertexID, initialGain);
        mutex.unlock();
    }

}

void updateBoundaryVertices(Graph& graph, std::vector<int>& boundaryVertices, std::vector<std::vector<int>>& initial_partitions, int node){
     std::vector<int> neighbors = graph.getNeighbors(node);

     std::vector<int> neighborNeigbhors;
     for(int vertex : neighbors){
         neighborNeigbhors = graph.getNeighbors(vertex);

         // Get the partition of the current vertex
         int partitionVertex = -1; // Initialize to an invalid partition

         // Find the partition that contains the current vertex
         for (int i = 0; i < initial_partitions.size(); ++i) {
             if (std::find(initial_partitions[i].begin(), initial_partitions[i].end(), vertex) != initial_partitions[i].end()) {
                 partitionVertex = i;
                 break;
             }
         }

         // Initialize a flag to check if the vertex is a boundary vertex
         bool isBoundary = false;

         // Iterate through the neighbors of the current vertex
         for (int neighbor : neighborNeigbhors) {
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

         bool wasBoundary = std::find(boundaryVertices.begin(), boundaryVertices.end(), vertex) != boundaryVertices.end();

         // If the vertex is a boundary vertex, add it to the list
         if (isBoundary) {
             if(!wasBoundary){
                 boundaryVertices.push_back(vertex);
             }
         }else{
            if(wasBoundary){
                boundaryVertices.erase(std::remove(boundaryVertices.begin(), boundaryVertices.end(), vertex), boundaryVertices.end());
            }
         }
     }

}

void balancePartitions(Graph& graph, std::vector<std::vector<int>>& partitions, std::vector<double>& partition_weights, std::vector<int>& boundaryVertices, std::vector<std::vector<int>> initial_partitions, float maxDeviation, double target_weight){
    int iterations = 0;
    //Perform movements until is balanced, or until it iterates one time for vertex
    while(iterations < MAX_SWAPS_PER_ITERATION){
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
                updateBoundaryVertices(graph, boundaryVertices, initial_partitions, node);
                break;
            }

            for(int vertex : partitions[partitionOrder[0]]){
                std::vector<int>& neighbors = graph.getNeighbors(vertex);

                for(int node : partitions[partitionOrder[i]]){
                    //Check all nodes of the heaviest partition not already checked
                    for(int neighbor : neighbors){
                        if(neighbor == node){
                            //move neighbor to partition
                            partitions[partitionOrder[i]].erase(std::remove(partitions[partitionOrder[i]].begin(), partitions[partitionOrder[i]].end(), neighbor), partitions[partitionOrder[i]].end());
                            partitions[partitionOrder[0]].push_back(node);
                            partition_weights[partitionOrder[i]] -= graph.getVertexWeight(neighbor);
                            partition_weights[partitionOrder[0]] += graph.getVertexWeight(neighbor);
                            swapped = true;
                            updateBoundaryVertices(graph, boundaryVertices, initial_partitions, neighbor);
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
}

std::vector<std::vector<int>> Uncoarsening(Graph& graph, std::vector<std::vector<int>>& partitions, std::vector<int>& boundaryVertices, int numThreads, float maxDeviation){

    std::cout << "Uncoarsening entered" << std::endl;
    int num_levels = graph.getCoarsingLevel();

    std::unordered_map<int, int> coarser_to_finer_mapping;
    std::unordered_map<std::pair<int, int>, double, HashPair> edgesMapping;
    std::unordered_map<int, double> verticesWeights;

    // Create a new vector to store the uncoarsed partitions
    std::vector<std::vector<int>> uncoarsed_partitions;

    // Create a new graph to store the uncoarsened version
    Graph uncoarsened_temp;

    //Create a new vector to store uncoarsened boundaries
    std::vector<int> uncoarsedBoundaries;

    //Partition weights
    std::vector<double> partition_weights(partitions.size(), 0.0);
    for (size_t i = 0; i < partitions.size(); ++i) {
        for (int element : partitions[i]) {
            partition_weights[i] += graph.getVertexWeight(element);
        }
    }

    double target_weight = 0.0;
    for(double weight : partition_weights){
        target_weight += weight;
    }
    target_weight = target_weight/partition_weights.size();

    //Data structure to contain vertex gains
    std::vector<std::pair<int, double>> vertexGains;

    for (int level = num_levels - 1; level >= 0; level--) {
        std::cout << "Uncoarsening a new level" << std::endl;
        // Create a single mutex for synchronization
        std::mutex mutex;

        //Get this level structures
        coarser_to_finer_mapping = graph.getMapping(level);
        edgesMapping = graph.getMappingEdges(level);
        verticesWeights = graph.getMappingVerticesWeights(level);

        //Graph uncoarsening
        std::cout << "Graph uncoarsening" << std::endl;

        // Divide the work among three threads
        int mappingPerThread = coarser_to_finer_mapping.size() / numThreads;
        int edgesPerThread = edgesMapping.size() / numThreads;

        //Create a barrier mutex
        std::mutex barrierMutex;
        crrnt_ctr = 0;  //reset barrier ctr
        std::condition_variable cv;

        // Create and join threads in a loop
        std::vector<std::thread> threads;
        for (int i = 0; i < numThreads; i++) {
            mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
            int start_mapping = i * mappingPerThread;
            int end_mapping = (i + 1) * mappingPerThread;
            int start_edges = i * edgesPerThread;
            int end_edges = (i + 1) * edgesPerThread;

            if (i == numThreads - 1) {
                // Handle the last thread's range
                end_mapping = coarser_to_finer_mapping.size();
                end_edges = edgesMapping.size();
            }

            threads.emplace_back([&] {
                UncoarseningGraph(uncoarsened_temp, coarser_to_finer_mapping, edgesMapping, verticesWeights, std::ref(mutex), start_mapping, end_mapping, start_edges, end_edges, numThreads, std::ref(cv), std::ref(barrierMutex));
            });
        }
        // Wait for all threads to finish
        for (std::thread& thread : threads) {
            thread.join();
        }
        threads.clear();
        uncoarsened_temp.copyCoarseningData(graph);
        graph = uncoarsened_temp;
        uncoarsened_temp.clear();
        std::cout << "Graph uncoarsening finished" << std::endl;

        //Uncoarse partitions
        std::cout << "Partition uncoarsening started" << std::endl;
        uncoarsed_partitions = std::vector<std::vector<int>>(partitions.size());

        // Divide the work among three threads
        int partitionsPerThread = coarser_to_finer_mapping.size() / numThreads;

        // Create and join threads in a loop
        for (int i = 0; i < numThreads; i++) {
            mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
            int start = i * partitionsPerThread;
            int end = (i + 1) * partitionsPerThread;

            if (i == numThreads - 1) {
                end = coarser_to_finer_mapping.size(); // Handle the last thread's range
            }

            threads.emplace_back([&] {
                UncoarsePartitions(coarser_to_finer_mapping, partitions, uncoarsed_partitions, std::ref(mutex), start, end);
            });
        }
        // Wait for all threads to finish
        for (std::thread& thread : threads) {
            thread.join();
        }
        partitions = uncoarsed_partitions;
        uncoarsed_partitions.clear();
        threads.clear();
        std::cout << "Uncoarse prtition finished" << std::endl;

        //Uncoarse boundary vertices

        std::cout << "Uncoarse boundary vertices started" << std::endl;
        //initialize uncoarsed boundary vertices
        uncoarsedBoundaries.clear();

        // Divide the work among three threads
        int verticesPerThread = boundaryVertices.size() / numThreads;

        // Create and join threads in a loop
        for (int i = 0; i < numThreads; i++) {
            mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
            int start = i * verticesPerThread;
            int end = (i + 1) * verticesPerThread;

            if (i == numThreads - 1) {
                end = graph.getVertices().size(); // Handle the last thread's range
            }

            threads.emplace_back([&] {
                uncoarseBoundaries(graph, coarser_to_finer_mapping, boundaryVertices, uncoarsedBoundaries, partitions, std::ref(mutex), start, end);
            });
        }
        // Wait for all threads to finish
        for (std::thread& thread : threads) {
            thread.join();
        }
        boundaryVertices = uncoarsedBoundaries;
        uncoarsedBoundaries.clear();
        threads.clear();

        std::cout << "Uncoarse boundary vertices finished" << std::endl;

        //If necessary, balance partitions
        std::cout << "Balance partitions" << std::endl;
        balancePartitions(graph, partitions, partition_weights, boundaryVertices, partitions, maxDeviation, target_weight);
        std::cout << "Balance partitions ended" << std::endl;

        //Refinement step
        std::cout << "Refinement step started" << std::endl;
        // Divide the work among three threads
        int vertexPerThread = boundaryVertices.size() / numThreads;

        // Create and join threads in a loop
        for (int i = 0; i < numThreads; i++) {
            mutex.lock();  //This will be unlocked by thread after it is correctly created, to avoid to pass to multiple threds the same start and end
            int start = i * vertexPerThread;
            int end = (i + 1) * vertexPerThread;

            if (i == numThreads - 1) {
                end = boundaryVertices.size(); // Handle the last thread's range
            }

            threads.emplace_back([&] {
                refinementStep(graph, vertexGains, partitions, boundaryVertices, std::ref(mutex), start, end);
            });
        }
        // Wait for all threads to finish
        for (std::thread& thread : threads) {
            thread.join();
        }
        boundaryVertices = uncoarsedBoundaries;
        uncoarsedBoundaries.clear();
        threads.clear();

        std::cout << "Vertex gain calculated" << std::endl;

        // Sort in descending order by gain vertexGain
        std::sort(vertexGains.begin(), vertexGains.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
        });

        std::vector<double> bestPartition;

        int swaps = 0;

        std::cout << "Start swapping" << std::endl;

        // Access pairs sequentially while gain > 0
        for (const auto& pair : vertexGains) {
            swaps++;
            if (pair.second <= 0.0 || swaps > MAX_SWAPS_PER_ITERATION) {
                break; // Stop when gain is not positive
            }

            // To reset the vector to the initial configuration with all elements set to -1:
            bestPartition = std::vector<double>(partitions.size(), 0.0);
            int vertexID = pair.first;

            // Get the partition of the current vertex
            int partitionVertex = -1; // Initialize to an invalid partition

            // Find the partition that contains the current vertex
            for (int i = 0; i < partitions.size(); ++i) {
                if (std::find(partitions[i].begin(), partitions[i].end(), vertexID) != partitions[i].end()) {
                    partitionVertex = i;
                    break;
                }
            }

            //Get neighbors
            std::vector<int>& neighbors = graph.getNeighbors(vertexID);

            for(int neighbor : neighbors){
                double weight = graph.getEdgeWeight(vertexID, neighbor);

                // Find the partition that contains the neighbor
                int partitionNeighbor = -1; // Initialize to an invalid partition
                for (int i = 0; i < partitions.size(); ++i) {
                    if (std::find(partitions[i].begin(), partitions[i].end(), neighbor) != partitions[i].end()) {
                        partitionNeighbor = i;
                        break;
                    }
                }
                bestPartition[partitionNeighbor] += weight;
            }

            //find the best partition
            std::vector<int> partitionOrder(bestPartition.size());

            // Initialize partitionOrder to contain indices from 0 to n-1
            for (int i = 0; i < partitionOrder.size(); ++i) {
                partitionOrder[i] = i;
            }

            // Sort in incresing order
            std::sort(partitionOrder.begin(), partitionOrder.end(), [&](int a, int b) {
                return bestPartition[a] > bestPartition[b];
            });

            for(int i=0; i<bestPartition.size(); i++){
                if(partitionOrder[i] == partitionVertex || partition_weights[partitionOrder[i]] == partition_weights[partitionVertex]){
                    break; //It is no more a gain
                }

                //Check if the swap would keep the balance
                double movedVertexWeight = graph.getVertexWeight(vertexID);
                double oldPartitionWeight = partition_weights[partitionVertex] - movedVertexWeight;
                double newPartitionWeight = partition_weights[partitionOrder[i]] + movedVertexWeight;

                if(oldPartitionWeight < maxDeviation*target_weight && oldPartitionWeight > (1 - (maxDeviation - 1))*target_weight && newPartitionWeight < maxDeviation*target_weight && newPartitionWeight > (1 - (maxDeviation - 1))*target_weight){
                    //move node to best partition
                    partitions[partitionVertex].erase(std::remove(partitions[partitionVertex].begin(), partitions[partitionVertex].end(), vertexID), partitions[partitionVertex].end());
                    partitions[partitionOrder[i]].push_back(vertexID);
                    partition_weights[partitionVertex] -= graph.getVertexWeight(vertexID);
                    partition_weights[partitionOrder[i]] += graph.getVertexWeight(vertexID);
                    updateBoundaryVertices(graph, boundaryVertices, partitions, vertexID);
                    break;
                }
            }
        }

        vertexGains.clear();
        std::cout << "Refinement step and cycle ended" << std::endl;
    }

    return partitions;
}