#ifndef MULKL_H
#define MULKL_H

// Enable or disable debug output
//#define DEBUG 1

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#include <iostream>
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

#include "KLCore.cpp"

// Structure to represent a partition
struct Partition {
    std::vector<int> part;
    double weight;
};

#ifdef DEBUG
// Function to calculate the cut size of partitions
double calculateCutSize(const Graph& G, const std::vector<std::vector<int>>& partitions);
// Function to print partitions
void printPartitions(const std::vector<std::vector<int>>& partitions, std::unordered_map<int, double> vertices);
#endif

// Function to fill a partition with dummy nodes to match the maximum weight
void fillPartition(Partition& partition, int& alreadyIns, const double& max_weight);
// Function to clean partitions by removing dummy nodes
std::vector<std::vector<int>> cleanPartition(const std::vector<Partition>& partitions);
// Function to create initial node partitions
std::vector<std::vector<int>> makeNodePartition_multiple_round_robin(const Graph& G, const int& numPartitions);
// Function to create initial node partitions
std::vector<std::vector<int>> makeNodePartition(Graph& G, bool expand);

std::pair<int, int> dividePartitionCount(const int &n);
// Function for Kernighan-Lin Multi-level Partitioning
std::vector<std::vector<int>> multipartitionKL_blob(Graph &G, int numPartitions);
// Function for Kernighan-Lin Multi-level Partitioning
std::vector<std::vector<int>> multipartitionKL_round_robin(Graph& G, const int& numPartitions);

#endif //MULKL_H