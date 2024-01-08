#define DEBUG 1

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

#include <chrono>
#include <iostream>
#include "ReaderWriter.cpp"
#include <random>



#define IN_FILE_NAME "resources\\metismodels\\x200000y440000m20q20.metis"
#define OUT_FILE_NAME "writingTest.prt"
#define NUM_THREADS 7

std::vector<std::vector<int>> generateRandomVector(const unsigned int& numPartitions, int N) {
    // Initialize a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    DEBUG_STDOUT("Number of Partitions: " + std::to_string(numPartitions));
    DEBUG_STDOUT("Number of Nodes: " + std::to_string(N));
    int nodesPerPartitions = N/numPartitions;
    DEBUG_STDOUT("Number of Nodes per Partition: " + std::to_string(nodesPerPartitions));
    // Random number of elements for each vector
    std::uniform_int_distribution<int> partitionChoice(0, numPartitions-1);

    // Generate the numbers
    std::vector<int> numbers;
    for(int i=1; i <= N; ++i){
        numbers.push_back(i);
    }
    // Vector of vectors to store the random numbers
    DEBUG_STDOUT("Size of numbers: "+ std::to_string(numbers.size()));

    std::vector<std::vector<int>> result;
    // initiallize inner Vectors
    for(unsigned int i = 0; i< numPartitions; i++)
        result.push_back({});

    int index;
    for (int n: numbers) {
        do{
            index = partitionChoice(gen);
        }while (result[index].size() > static_cast<unsigned int> (nodesPerPartitions*3/2));

        result[index].push_back(n);
    }

    return result;
}



int main() {
    // Start the clock
    DEBUG_STDOUT("Started Reading");
    auto start_time = std::chrono::high_resolution_clock::now();
    
    //int nthreads = std::stoi(argv[1]);
    auto graph = metisRead(IN_FILE_NAME, NUM_THREADS);
    // Stop the clock
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    double seconds = duration.count() / 1e6;
    // Print the Reading execution time
    std::cout << "Reading execution time: " << seconds << " seconds" << std::endl;
    DEBUG_STDOUT("Nodes: " + std::to_string(graph.numVertices()) + "  Edges: " + std::to_string(graph.numEdges()));
    //graph.print();


    auto partitions = generateRandomVector(7, graph.numVertices());
    
    // // Display the generated vectors
    // for (const auto& part : partitions) {
    //     std::cout << "Vector: ";
    //     for (const int& node : part) {
    //         std::cout << node << " ";
    //     }
    //     std::cout << std::endl;
    // }


    DEBUG_STDOUT("Started Writing");
    // Start the clock
    start_time = std::chrono::high_resolution_clock::now();
    writeToFile(partitions, OUT_FILE_NAME);
    // Stop the clock
    end_time = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Convert the duration to a double value in seconds
    seconds = duration.count() / 1e6;
    // Print the Reading execution time
    std::cout << "Reading execution time: " << seconds << " seconds" << std::endl;
    
    return 0;
}

